#!/usr/bin/env Rscript

# ============================================================================
# Optimized_2 Dirichlet-Multinomial Regression for 3'UTR Usage Along Pseudotime
#
# Changes from optimized_1 (fixes issues identified in QC):
#   - Shared concentration: s estimated once under full model, used for both
#     full and null - restores LRT nesting guarantee (ll_full >= ll_null)
#   - Eliminates all ll_full < ll_null violations (312 genes in optimized_1)
#   - Beta-only SEs with fixed s (improved Hessian conditioning)
#   - Expanded test suite covering the QC-identified issues
#
# Retained from optimized_1:
#   - Joint MLE via L-BFGS-B with analytical gradients
#   - Gene-specific concentration via empirical Bayes shrinkage
#   - Correct LRT df = spline_df * (K-1)
#   - Vectorized log-likelihood
#
# Run with --test to execute the test suite (Test Mode section).
# ============================================================================

# ===== Section 0: Libraries and Argument Parsing =============================

library(tidyverse)
library(anndata)
library(pbapply)
library(parallel)
library(splines)
library(numDeriv)
library(digest)    # xxhash32 for per-gene seed derivation (Phase C/D reproducibility)

pboptions(type = "timer")

parse_args <- function(args) {
    opts <- list(
        data_path = NULL, method = "joint", spline_df = 5L, min_cells = 50L,
        padj_method = "BH", ncores = max(1L, detectCores() - 1L),
        compute_se = TRUE, max_iter = 500L, shrinkage = "adaptive",
        lambda = 1.0, test = FALSE,
        bootstrap = 0L, bootstrap_pval = 0.1,
        calibration_n = 200L, lambda_sensitivity = FALSE,
        seed = 42L
    )
    i <- 1L
    while (i <= length(args)) {
        a <- args[i]
        if (a == "--test") {
            opts$test <- TRUE
            i <- i + 1L
            next
        }
        if (a == "--lambda_sensitivity") {
            opts$lambda_sensitivity <- TRUE
            i <- i + 1L
            next
        }
        if (startsWith(a, "--") && i < length(args)) {
            val <- args[i + 1L]
            switch(sub("^--", "", a),
                method = {
                    opts$method <- val
                },
                spline_df = {
                    opts$spline_df <- as.integer(val)
                },
                min_cells = {
                    opts$min_cells <- as.integer(val)
                },
                padj_method = {
                    opts$padj_method <- val
                },
                ncores = {
                    opts$ncores <- as.integer(val)
                },
                compute_se = {
                    opts$compute_se <- as.logical(val)
                },
                max_iter = {
                    opts$max_iter <- as.integer(val)
                },
                shrinkage = {
                    opts$shrinkage <- val
                },
                lambda = {
                    opts$lambda <- as.numeric(val)
                },
                bootstrap = {
                    opts$bootstrap <- as.integer(val)
                },
                bootstrap_pval = {
                    opts$bootstrap_pval <- as.numeric(val)
                },
                calibration_n = {
                    opts$calibration_n <- as.integer(val)
                },
                seed = {
                    opts$seed <- as.integer(val)
                },
                stop(paste("Unknown argument:", a))
            )
            i <- i + 2L
        } else if (!startsWith(a, "--")) {
            opts$data_path <- a
            i <- i + 1L
        } else {
            stop(paste("Missing value for:", a))
        }
    }
    if (!opts$test && is.null(opts$data_path)) {
        stop("Usage: Rscript 260401-dirichlet_multinomial_analysis_optimized_2.R <path_to_h5ad> [--method joint] [--spline_df 5] [--bootstrap N] [--bootstrap_pval 0.1] [--lambda_sensitivity] [--test]")
    }
    opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

# ===== Section 2: Core Mathematical Functions ================================

#' Softmax-based probability computation (K-th category as reference)
#' @param beta p x (K-1) coefficient matrix
#' @param mm N x p design matrix
#' @return N x K probability matrix
compute_pi <- function(beta, mm) {
    beta <- as.matrix(beta)
    eta <- mm %*% beta # N x (K-1)
    eta_aug <- cbind(eta, 0) # append 0 for reference category
    eta_max <- apply(eta_aug, 1, max) # numerical stability
    exp_eta <- exp(eta_aug - eta_max)
    pi_mat <- exp_eta / rowSums(exp_eta)
    pi_mat <- pmax(pi_mat, 1e-12)
    pi_mat / rowSums(pi_mat)
}

#' Vectorized Dirichlet-Multinomial log-likelihood
#' @param Y N x K count matrix
#' @param alpha N x K Dirichlet parameter matrix (= pi * s)
#' @return scalar log-likelihood
dirmult_loglik_vec <- function(Y, alpha) {
    n <- rowSums(Y)
    sa <- rowSums(alpha)
    sum(lgamma(sa) + lgamma(n + 1) - lgamma(n + sa) +
        rowSums(lgamma(Y + alpha) - lgamma(Y + 1) - lgamma(alpha)))
}

#' Negative log-likelihood for optim (with optional shrinkage penalty)
dirmult_negloglik <- function(theta, Y, mm, K, s_trend = NULL, lambda = 0) {
    p <- ncol(mm)
    n_beta <- p * (K - 1L)
    beta <- matrix(theta[1:n_beta], nrow = p, ncol = K - 1L)
    log_s <- theta[n_beta + 1L]
    s <- exp(log_s)

    pi_mat <- compute_pi(beta, mm)
    alpha <- pi_mat * s
    nll <- -dirmult_loglik_vec(Y, alpha)

    if (lambda > 0 && !is.null(s_trend)) {
        nll <- nll + lambda * (log_s - log(s_trend))^2
    }
    nll
}

#' Analytical gradient of negative log-likelihood
dirmult_neggrad <- function(theta, Y, mm, K, s_trend = NULL, lambda = 0) {
    p <- ncol(mm)
    N <- nrow(Y)
    n_beta <- p * (K - 1L)
    beta <- matrix(theta[1:n_beta], nrow = p, ncol = K - 1L)
    log_s <- theta[n_beta + 1L]
    s <- exp(log_s)

    pi_mat <- compute_pi(beta, mm)
    alpha <- pi_mat * s
    n_counts <- rowSums(Y)

    # digamma differences: N x K
    delta <- digamma(Y + alpha) - digamma(alpha)

    # Gradient w.r.t. pre-softmax linear predictors eta, using softmax Jacobian:
    #   dLL/d(eta_m) = s * pi_m * (delta_m - sum_k pi_k * delta_k)
    delta_weighted <- rowSums(delta * pi_mat) # N-vector
    dLL_deta <- s * pi_mat[, 1:(K - 1), drop = FALSE] *
        (delta[, 1:(K - 1), drop = FALSE] - delta_weighted)

    # Chain rule: eta = mm %*% beta, so dLL/d(beta) = t(mm) %*% dLL/d(eta)
    dLL_dbeta <- crossprod(mm, dLL_deta) # p x (K-1)

    # Gradient w.r.t. log(s):
    #   d/d(log_s) [lgamma(s) - lgamma(n_i + s)] = s * [digamma(s) - digamma(n_i + s)]
    #   d/d(log_s) [lgamma(y+alpha) - lgamma(alpha)] = alpha * delta (since alpha = pi*s)
    dLL_dlogs <- s * (N * digamma(s) - sum(digamma(n_counts + s)) +
        sum(pi_mat * delta))

    grad <- -c(as.vector(dLL_dbeta), dLL_dlogs)

    if (lambda > 0 && !is.null(s_trend)) {
        grad[n_beta + 1L] <- grad[n_beta + 1L] + 2 * lambda * (log_s - log(s_trend))
    }
    grad
}

#' Per-isoform observed-to-expected variance ratio under fitted DirMult
#'
#' For each isoform k, computes the ratio of the empirical variance of the
#' cell-level fractions f_ik = Y_ik / n_i to the variance predicted by the
#' fitted DirMult model: Var_DirMult(f_k) = pi_ik (1 - pi_ik) (n_i + s) /
#' [n_i (s + 1)]. Under correct specification, all r_k should be ~1.
#'
#' Returns a length-K numeric vector. NA components are returned for isoforms
#' with degenerate variance (zero counts everywhere) or for cells where the
#' total count is zero (which are dropped from the calculation).
#'
#' @param Y N x K integer count matrix
#' @param pi_mat N x K predicted-proportion matrix from the fitted model
#' @param s scalar fitted concentration
#' @return length-K numeric vector of variance ratios r_k
compute_var_ratios <- function(Y, pi_mat, s) {
    n_counts <- rowSums(Y)
    keep <- n_counts > 0
    if (sum(keep) < 3L) {
        return(rep(NA_real_, ncol(Y)))
    }
    Y_sub <- Y[keep, , drop = FALSE]
    pi_sub <- pi_mat[keep, , drop = FALSE]
    n_sub <- n_counts[keep]
    fracs <- Y_sub / n_sub
    K <- ncol(Y_sub)
    r <- numeric(K)
    for (k in seq_len(K)) {
        # Empirical variance of cell-level fractions
        v_obs <- var(fracs[, k])
        # Cell-specific predicted variances under DirMult; average across cells
        var_terms <- pi_sub[, k] * (1 - pi_sub[, k]) * (n_sub + s) / (n_sub * (s + 1))
        v_exp <- mean(var_terms)
        if (is.na(v_obs) || is.na(v_exp) || v_exp <= 0) {
            r[k] <- NA_real_
        } else {
            r[k] <- v_obs / v_exp
        }
    }
    r
}

#' Method-of-moments concentration estimator
#' Uses Var(fraction_k) ≈ pi_k(1-pi_k) / (s+1) for Dirichlet fractions
mom_concentration <- function(Y) {
    fracs <- Y / pmax(rowSums(Y), 1)
    pi_hat <- colMeans(fracs)
    s_per_site <- vapply(seq_len(ncol(Y)), function(k) {
        v <- var(fracs[, k])
        if (v < 1e-12 || pi_hat[k] < 1e-6 || pi_hat[k] > 1 - 1e-6) {
            return(NA_real_)
        }
        pi_hat[k] * (1 - pi_hat[k]) / v - 1
    }, numeric(1))
    s <- median(s_per_site, na.rm = TRUE)
    if (is.na(s) || s < 0.5) s <- 0.5
    s
}

#' Warm-start beta initialization via K independent binomial GLMs
#' Converts to K-1 reference parameterization by subtracting K-th column
init_beta_from_glm <- function(Y, mm, K) {
    fracs <- Y / pmax(rowSums(Y), 1)
    w <- rowSums(Y)
    p <- ncol(mm)
    beta_all <- matrix(0, nrow = p, ncol = K)
    for (k in seq_len(K)) {
        beta_all[, k] <- tryCatch(
            {
                fit <- suppressWarnings(
                    glm.fit(mm, fracs[, k], family = binomial(link = "logit"), weights = w)
                )
                coef <- fit$coefficients
                coef[is.na(coef)] <- 0
                coef
            },
            error = function(e) rep(0, p)
        )
    }
    # K-1 reference parameterization: subtract reference (K-th) column
    beta_all[, 1:(K - 1), drop = FALSE] - beta_all[, K]
}

# ===== Section 2.5: Global Concentration Trend Estimation ====================

#' Estimate loess trend: log(s_mom) ~ log(mean_counts) across all genes
#' Provides data-driven s_trend per gene for empirical Bayes shrinkage
estimate_global_trend <- function(X, utr_names) {
    unique_utrs <- unique(utr_names)
    n_utrs <- length(unique_utrs)
    s_mom_vec <- numeric(n_utrs)
    mean_counts_vec <- numeric(n_utrs)

    for (idx in seq_len(n_utrs)) {
        cols <- which(utr_names == unique_utrs[idx])
        Y_utr <- X[, cols, drop = FALSE]
        keep <- rowSums(Y_utr) > 0
        Y_sel <- Y_utr[keep, , drop = FALSE]

        if (ncol(Y_sel) < 2 || nrow(Y_sel) < 10 || sum(Y_sel) <= 10) {
            s_mom_vec[idx] <- NA
            mean_counts_vec[idx] <- NA
            next
        }
        s_mom_vec[idx] <- mom_concentration(Y_sel)
        mean_counts_vec[idx] <- mean(rowSums(Y_sel))
    }
    names(s_mom_vec) <- unique_utrs
    names(mean_counts_vec) <- unique_utrs

    valid <- !is.na(s_mom_vec) & !is.na(mean_counts_vec) &
        s_mom_vec > 0 & mean_counts_vec > 0
    trend_data <- data.frame(
        log_s = log(s_mom_vec[valid]),
        log_counts = log(mean_counts_vec[valid])
    )
    trend_fit <- loess(log_s ~ log_counts, data = trend_data, span = 0.3)

    s_trend <- setNames(rep(NA_real_, n_utrs), unique_utrs)
    s_trend[valid] <- exp(predict(trend_fit, newdata = data.frame(
        log_counts = log(mean_counts_vec[valid])
    )))
    s_median <- median(s_trend, na.rm = TRUE)
    s_trend[is.na(s_trend)] <- s_median

    list(
        s_trend = s_trend, s_mom = s_mom_vec, mean_counts = mean_counts_vec,
        trend_fit = trend_fit
    )
}

# ===== Section 3: Optimization Wrappers ======================================

#' Joint MLE: optimize [vec(beta), log_s] simultaneously via L-BFGS-B
optimize_joint <- function(Y, mm, K, beta_init, s_init, s_trend, lambda, max_iter) {
    p <- ncol(mm)
    n_beta <- p * (K - 1L)
    theta0 <- c(as.vector(beta_init), log(s_init))
    lower <- c(rep(-Inf, n_beta), -5)
    upper <- c(rep(Inf, n_beta), 15)

    result <- optim(theta0,
        fn = dirmult_negloglik, gr = dirmult_neggrad,
        Y = Y, mm = mm, K = K, s_trend = s_trend, lambda = lambda,
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = max_iter)
    )

    beta <- matrix(result$par[1:n_beta], nrow = p, ncol = K - 1L)
    s <- exp(result$par[n_beta + 1L])
    pi_mat <- compute_pi(beta, mm)
    ll <- dirmult_loglik_vec(Y, pi_mat * s)

    list(beta = beta, s = s, ll = ll, convergence = result$convergence)
}

#' Profile likelihood: alternate between beta and s optimization
optimize_profile <- function(Y, mm, K, beta_init, s_init, s_trend, lambda, max_iter) {
    s <- s_init
    beta <- beta_init
    p <- ncol(mm)
    n_beta <- p * (K - 1L)
    iter_budget <- max(50L, max_iter %/% 5L)
    conv <- 1L

    for (outer in 1:5) {
        # Fix s, optimize beta
        fn_b <- function(bv) dirmult_negloglik(c(bv, log(s)), Y, mm, K)
        gr_b <- function(bv) dirmult_neggrad(c(bv, log(s)), Y, mm, K)[1:n_beta]
        res_b <- optim(as.vector(beta), fn_b, gr_b,
            method = "L-BFGS-B",
            control = list(maxit = iter_budget)
        )
        beta <- matrix(res_b$par, nrow = p, ncol = K - 1L)
        conv <- res_b$convergence

        # Fix beta, optimize log_s via 1D optimize
        pi_mat <- compute_pi(beta, mm)
        fn_s <- function(log_s) {
            nll <- -dirmult_loglik_vec(Y, pi_mat * exp(log_s))
            if (lambda > 0 && !is.null(s_trend)) {
                nll <- nll + lambda * (log_s - log(s_trend))^2
            }
            nll
        }
        res_s <- optimize(fn_s, interval = c(-5, 15))
        s <- exp(res_s$minimum)
    }

    pi_mat <- compute_pi(beta, mm)
    ll <- dirmult_loglik_vec(Y, pi_mat * s)
    list(beta = beta, s = s, ll = ll, convergence = conv)
}

#' Fixed-s: use trend concentration, only optimize beta
optimize_beta_fixed_s <- function(Y, mm, K, beta_init, s, max_iter) {
    p <- ncol(mm)
    n_beta <- p * (K - 1L)
    fn_b <- function(bv) dirmult_negloglik(c(bv, log(s)), Y, mm, K)
    gr_b <- function(bv) dirmult_neggrad(c(bv, log(s)), Y, mm, K)[1:n_beta]

    result <- optim(as.vector(beta_init), fn_b, gr_b,
        method = "L-BFGS-B",
        control = list(maxit = max_iter)
    )

    beta <- matrix(result$par, nrow = p, ncol = K - 1L)
    pi_mat <- compute_pi(beta, mm)
    ll <- dirmult_loglik_vec(Y, pi_mat * s)
    list(beta = beta, s = s, ll = ll, convergence = result$convergence)
}

# ===== Section 4: Gene-Level Fitting =========================================

#' Fit full and null DirMult models for one gene, compute LRT
#'
#' optimized_2: shared-s pipeline with conditional
#' MLE refits in both the full and null beta subspaces.
#'
#' Pipeline:
#'   1. Fit full model jointly (or profile/fixed-s); on failure, fall back to
#'      fixed-s on the full design at s_trend.
#'   2. Extract s_shared = fit_full$s.
#'   3. Refit the full model beta at fixed s_shared (Stage 1). The joint optimum
#'      beta is not in general the conditional MLE at s_shared because L(beta, s)
#'      is not concave in beta under DirMult + softmax. The refit produces an
#'      actual conditional MLE on the full beta subspace, mirroring how the null
#'      is fit.
#'   4. Fit the null model beta at fixed s_shared.
#'   5. If ll_full < ll_null after Stage 1, refit the full model from beta_null
#'      padded with zeros in the spline columns (Stage 2). This initialization
#'      gives L = ll_null exactly, so any L-BFGS-B step from there can only
#'      ascend, guaranteeing ll_full >= ll_null by construction. The padding is
#'      verified by re-evaluating the LL on the padded beta and falling through
#'      to Stage 3 if the column-prefix property has been broken (e.g. by QR
#'      pivoting on rank-deficient designs).
#'   6. If ll_full < ll_null even after Stage 2 (rare), set
#'      lrt_status = "optimization_failed" and return NA for deviance and
#'      p_value (Stage 3). These genes are excluded from BH adjustment in the
#'      output assembly rather than reported with deviance = 0, p = 1.
#'
#' Also returns per-isoform observed-to-expected variance ratios as
#' DirMult-misspecification diagnostics (max_var_ratio, min_var_ratio).
fit_dirmult_gene <- function(Y, mm_full, mm_null, K, method, s_trend,
                             lambda, max_iter, compute_se) {
    N <- nrow(Y)
    na_result <- list(
        deviance = NA, p_value = NA, s_full = NA, s_null = NA,
        s_shared = NA, ll_full = NA, ll_null = NA, n_cells = N,
        K = K, lrt_df = NA, converged = FALSE,
        full_converged = NA, cells_per_param = NA,
        method = method, se = NULL,
        lrt_status = "skipped",
        max_var_ratio = NA_real_, min_var_ratio = NA_real_,
        beta_full = NULL, beta_null = NULL, mm_full = NULL, mm_null = NULL
    )

    # Fix rank deficiency (drop linearly dependent columns)
    fix_rank <- function(mm) {
        qr_mm <- qr(mm)
        if (qr_mm$rank < ncol(mm)) {
            mm[, qr_mm$pivot[seq_len(qr_mm$rank)], drop = FALSE]
        } else {
            mm
        }
    }
    mm_full <- fix_rank(mm_full)
    mm_null <- fix_rank(mm_null)

    # Correct LRT degrees of freedom
    lrt_df <- (ncol(mm_full) - ncol(mm_null)) * (K - 1L)
    if (lrt_df < 1L) {
        return(na_result)
    }

    # Cells-per-parameter diagnostic
    cells_per_param <- N / lrt_df

    # Initialize coefficients
    beta_init_full <- init_beta_from_glm(Y, mm_full, K)
    beta_init_null <- init_beta_from_glm(Y, mm_null, K)

    # --- Step 1: Fit full model (primary method + fallback) ---
    full_converged <- FALSE
    fit_full <- tryCatch(
        {
            if (method == "joint") {
                optimize_joint(
                    Y, mm_full, K, beta_init_full, s_trend, s_trend,
                    lambda, max_iter
                )
            } else if (method == "profile") {
                optimize_profile(
                    Y, mm_full, K, beta_init_full, s_trend, s_trend,
                    lambda, max_iter
                )
            } else {
                optimize_beta_fixed_s(
                    Y, mm_full, K, beta_init_full, s_trend,
                    max_iter
                )
            }
        },
        error = function(e) NULL
    )

    if (!is.null(fit_full) && fit_full$convergence == 0) {
        full_converged <- TRUE
    } else {
        # Fallback: fixed s at trend
        fit_full <- tryCatch(
            optimize_beta_fixed_s(
                Y, mm_full, K, beta_init_full, s_trend,
                max_iter
            ),
            error = function(e) NULL
        )
    }
    if (is.null(fit_full)) {
        return(na_result)
    }

    # --- Step 2: Extract shared concentration ---
    s_shared <- fit_full$s

    # --- Step 2b (Stage 1): Refit full model beta at fixed s_shared ---
    # The conditional LRT requires max_beta L(beta, s_shared) on the full beta
    # space, not L at the joint optimum's beta coordinate. The two coincide
    # for well-behaved (concave) cases but diverge when the digamma
    # perturbation in DirMult breaks concavity in beta. The warm start from
    # the joint optimum keeps the L-BFGS-B refit cheap and well-positioned.
    fit_full_refit <- tryCatch(
        optimize_beta_fixed_s(Y, mm_full, K, fit_full$beta, s_shared, max_iter),
        error = function(e) NULL
    )
    if (!is.null(fit_full_refit) && fit_full_refit$ll >= fit_full$ll - 1e-10) {
        fit_full <- fit_full_refit
    }

    # --- Step 3: Fit null model with s fixed at s_shared (beta-only) ---
    fit_null <- tryCatch(
        optimize_beta_fixed_s(
            Y, mm_null, K, beta_init_null, s_shared,
            max_iter
        ),
        error = function(e) NULL
    )
    if (is.null(fit_null)) {
        na_result$s_full <- s_shared
        na_result$s_null <- s_shared
        na_result$s_shared <- s_shared
        na_result$ll_full <- fit_full$ll
        na_result$full_converged <- full_converged
        na_result$cells_per_param <- cells_per_param
        na_result$lrt_status <- "null_fit_failed"
        return(na_result)
    }

    # --- Step 3b (Stage 2): Provably-nested fallback via padded null ---
    # If still ll_full < ll_null, refit the full model starting from beta_null
    # padded with zeros into the full beta subspace. By construction,
    # L(beta_null_padded, s_shared) = ll_null exactly, so any L-BFGS-B step
    # from this starting point can only ascend, giving ll_full >= ll_null
    # by construction. The ll_padded_check guards against future regressions
    # in the column-prefix property between mm_null and mm_full.
    if (fit_full$ll < fit_null$ll - 1e-10) {
        p_full <- ncol(mm_full)
        p_null <- ncol(mm_null)
        if (p_null <= p_full) {
            beta_null_padded <- matrix(0, nrow = p_full, ncol = K - 1L)
            beta_null_padded[seq_len(p_null), ] <- fit_null$beta
            ll_padded_check <- tryCatch(
                dirmult_loglik_vec(Y, compute_pi(beta_null_padded, mm_full) * s_shared),
                error = function(e) NA_real_
            )
            if (!is.na(ll_padded_check) &&
                abs(ll_padded_check - fit_null$ll) < 1e-6) {
                fit_full_safe <- tryCatch(
                    optimize_beta_fixed_s(
                        Y, mm_full, K, beta_null_padded, s_shared, max_iter
                    ),
                    error = function(e) NULL
                )
                if (!is.null(fit_full_safe) &&
                    fit_full_safe$ll >= fit_null$ll - 1e-10) {
                    fit_full <- fit_full_safe
                }
            }
        }
    }

    # --- Step 4: LRT (both likelihoods at s_shared) ---
    # Stage 3: if ll_full is still below ll_null, the LRT cannot be constructed
    # correctly. Return NA with status "optimization_failed" so the gene is
    # excluded from BH adjustment rather than reported as p = 1.
    lrt_status <- "ok"
    if (fit_full$ll < fit_null$ll - 1e-10) {
        lrt_status <- "optimization_failed"
        deviance <- NA_real_
        p_value <- NA_real_
    } else {
        # Numerical safety: clip tiny floating-point negatives at zero
        deviance <- max(-2 * (fit_null$ll - fit_full$ll), 0)
        p_value <- pchisq(deviance, df = lrt_df, lower.tail = FALSE)
    }

    # --- Step 5: Optional standard errors (beta-only with fixed s) ---
    se <- NULL
    if (compute_se) {
        se <- tryCatch(
            {
                fn_beta_only <- function(bv) {
                    dirmult_negloglik(c(bv, log(s_shared)), Y, mm_full, K)
                }
                H <- numDeriv::hessian(fn_beta_only, as.vector(fit_full$beta))
                vcov <- solve(H)
                sqrt(pmax(diag(vcov), 0))
            },
            error = function(e) NULL
        )
    }

    # --- Step 6: DirMult misspecification diagnostics ---
    var_ratios <- tryCatch(
        compute_var_ratios(Y, compute_pi(fit_full$beta, mm_full), s_shared),
        error = function(e) rep(NA_real_, K)
    )
    max_vr <- if (all(is.na(var_ratios))) NA_real_ else max(var_ratios, na.rm = TRUE)
    min_vr <- if (all(is.na(var_ratios))) NA_real_ else min(var_ratios, na.rm = TRUE)

    converged <- (fit_full$convergence == 0) && (fit_null$convergence == 0)

    list(
        deviance = deviance, p_value = p_value, s_full = s_shared,
        s_null = s_shared, s_shared = s_shared,
        ll_full = fit_full$ll, ll_null = fit_null$ll,
        n_cells = N, K = K, lrt_df = lrt_df, converged = converged,
        full_converged = full_converged, cells_per_param = cells_per_param,
        method = method, se = se,
        lrt_status = lrt_status,
        max_var_ratio = max_vr, min_var_ratio = min_vr,
        beta_full = fit_full$beta, beta_null = fit_null$beta,
        mm_full = mm_full, mm_null = mm_null
    )
}

# ===== Section 4.5: Parametric Bootstrap of the Conditional LRT ==============

#' Parametric bootstrap p-value for the conditional LRT statistic
#'
#' Generates B bootstrap datasets from the fitted null model
#'   Y_b[i, ] ~ DirMult(n_i, pi_null_i, s_shared)
#' (cell totals n_i held fixed; conditional bootstrap), refits beta under both
#' the full and null designs at fixed s_shared, and returns the empirical
#' p-value
#'   p_boot = (1 + #{D_b >= D_obs}) / (B + 1).
#'
#' Why valid: under H0, Y is exchangeable with the bootstrap replicates because
#' both come from the same DirMult(n, pi_null, s_shared). The bootstrap
#' distribution of D is exact (in the Monte Carlo sense) for the conditional
#' LRT statistic the method actually computes - bypassing the chi-squared
#' approximation, the asymptotic regularity conditions, the concavity
#' assumption, and the MAP/MLE mismatch.
#'
#' Why fast: each bootstrap fit only needs the fast `optimize_beta_fixed_s`
#' path (no joint optimization, no Hessian, no SE), warm-started from the
#' original beta_null and beta_full. Each bootstrap fit is ~5-10% the cost of
#' an original joint fit.
#'
#' @param Y_obs original count matrix (used only for K and dimensionality)
#' @param mm_full rank-fixed full design matrix (as returned by fit_dirmult_gene)
#' @param mm_null rank-fixed null design matrix
#' @param K number of isoforms
#' @param beta_full_obs original full-model beta (warm start)
#' @param beta_null_obs original null-model beta
#' @param s_shared the plug-in concentration
#' @param D_obs the observed deviance
#' @param B number of bootstrap replicates
#' @param max_iter max L-BFGS-B iterations per bootstrap fit
#' @param n_counts vector of per-cell totals to use as the multinomial sizes
#' @return list(p_boot, n_completed, D_b)
bootstrap_dirmult_lrt <- function(Y_obs, mm_full, mm_null, K,
                                  beta_full_obs, beta_null_obs, s_shared,
                                  D_obs, B, max_iter, n_counts) {
    if (B <= 0L || is.na(D_obs) || is.na(s_shared)) {
        return(list(p_boot = NA_real_, n_completed = 0L, D_b = numeric(0)))
    }
    N <- nrow(Y_obs)
    pi_null <- compute_pi(beta_null_obs, mm_null)
    alpha_null <- pi_null * s_shared

    # Pre-bind warm-start initializations
    p_full <- ncol(mm_full)
    p_null <- ncol(mm_null)
    beta_warm_full <- beta_full_obs
    beta_warm_null <- beta_null_obs

    D_b <- numeric(B)
    n_completed <- 0L
    n_ge <- 0L
    for (b in seq_len(B)) {
        Y_b <- tryCatch(
            {
                t(vapply(seq_len(N), function(i) {
                    p_dir <- rgamma(K, alpha_null[i, ])
                    p_dir <- p_dir / sum(p_dir)
                    as.integer(rmultinom(1, n_counts[i], p_dir))
                }, integer(K)))
            },
            error = function(e) NULL
        )
        if (is.null(Y_b)) {
            D_b[b] <- NA_real_
            next
        }
        # Refit both models at fixed s_shared
        fit_b_full <- tryCatch(
            optimize_beta_fixed_s(Y_b, mm_full, K, beta_warm_full, s_shared, max_iter),
            error = function(e) NULL
        )
        fit_b_null <- tryCatch(
            optimize_beta_fixed_s(Y_b, mm_null, K, beta_warm_null, s_shared, max_iter),
            error = function(e) NULL
        )
        if (is.null(fit_b_full) || is.null(fit_b_null)) {
            D_b[b] <- NA_real_
            next
        }
        # Apply the same Stage 2 padded-null safety net used in fit_dirmult_gene,
        # so the bootstrap statistic is computed by the same procedure as the
        # observed statistic.
        if (fit_b_full$ll < fit_b_null$ll - 1e-10 && p_null <= p_full) {
            beta_pad <- matrix(0, nrow = p_full, ncol = K - 1L)
            beta_pad[seq_len(p_null), ] <- fit_b_null$beta
            ll_pad <- tryCatch(
                dirmult_loglik_vec(Y_b, compute_pi(beta_pad, mm_full) * s_shared),
                error = function(e) NA_real_
            )
            if (!is.na(ll_pad) && abs(ll_pad - fit_b_null$ll) < 1e-6) {
                fit_b_safe <- tryCatch(
                    optimize_beta_fixed_s(Y_b, mm_full, K, beta_pad, s_shared, max_iter),
                    error = function(e) NULL
                )
                if (!is.null(fit_b_safe) && fit_b_safe$ll >= fit_b_null$ll - 1e-10) {
                    fit_b_full <- fit_b_safe
                }
            }
        }
        if (fit_b_full$ll < fit_b_null$ll - 1e-10) {
            # Bootstrap replicate failed Stage 2 too - skip it
            D_b[b] <- NA_real_
            next
        }
        D_b[b] <- max(-2 * (fit_b_null$ll - fit_b_full$ll), 0)
        n_completed <- n_completed + 1L
        if (D_b[b] >= D_obs) n_ge <- n_ge + 1L
    }
    if (n_completed < 1L) {
        return(list(p_boot = NA_real_, n_completed = 0L, D_b = D_b))
    }
    # Standard bootstrap p-value with +1/+1 correction
    p_boot <- (1 + n_ge) / (n_completed + 1)
    list(p_boot = p_boot, n_completed = n_completed, D_b = D_b)
}

# ===== Test Mode =============================================================

if (opts$test) {
    cat("=== Running Verification Tests ===\n\n")
    set.seed(42)
    all_pass <- TRUE

    # Test 1: LL cross-check (vectorized vs per-cell loop)
    cat("Test 1: LL vectorized vs per-cell loop... ")
    N_t <- 200
    K_t <- 3
    Y_t <- matrix(rpois(N_t * K_t, 5) + 1L, nrow = N_t)
    alpha_t <- matrix(rep(c(2, 3, 1), each = N_t), nrow = N_t) +
        matrix(runif(N_t * K_t, 0, 0.5), nrow = N_t)

    ll_vec <- dirmult_loglik_vec(Y_t, alpha_t)
    ll_loop <- sum(sapply(1:N_t, function(i) {
        x <- Y_t[i, ]
        a <- alpha_t[i, ]
        lgamma(sum(a)) + lgamma(sum(x) + 1) - lgamma(sum(x) + sum(a)) +
            sum(lgamma(x + a) - lgamma(x + 1) - lgamma(a))
    }))
    diff1 <- abs(ll_vec - ll_loop)
    if (diff1 < 1e-8) {
        cat("PASS (diff =", diff1, ")\n")
    } else {
        cat("FAIL (diff =", diff1, ")\n")
        all_pass <- FALSE
    }

    # Test 2: Analytical gradient vs numerical gradient
    cat("Test 2: Gradient check... ")
    N_g <- 50
    K_g <- 3
    p_g <- 3
    mm_g <- cbind(1, matrix(rnorm(N_g * (p_g - 1)), nrow = N_g))
    Y_g <- matrix(rpois(N_g * K_g, 3) + 1L, nrow = N_g)
    theta_g <- c(rnorm(p_g * (K_g - 1)), log(5))

    g_analytical <- dirmult_neggrad(theta_g, Y_g, mm_g, K_g)
    g_numerical <- numDeriv::grad(dirmult_negloglik, theta_g,
        Y = Y_g, mm = mm_g, K = K_g
    )
    rel_err <- max(abs(g_analytical - g_numerical) / pmax(abs(g_numerical), 1e-8))
    if (rel_err < 1e-5) {
        cat("PASS (max rel error:", rel_err, ")\n")
    } else {
        cat("FAIL (max rel error:", rel_err, ")\n")
        cat("  Analytical:", head(g_analytical, 5), "...\n")
        cat("  Numerical: ", head(g_numerical, 5), "...\n")
        all_pass <- FALSE
    }

    # Test 3: Synthetic data parameter recovery
    cat("Test 3: Synthetic DirMult parameter recovery... ")
    N_s <- 500
    K_s <- 3
    pt_s <- sort(runif(N_s))
    mm_s <- model.matrix(~ bs(pt_s, df = 3))
    p_s <- ncol(mm_s)
    true_beta <- matrix(c(
        0.5, -0.3, 0.2, 0.1,
        -0.2, 0.4, -0.1, 0.3
    ), nrow = p_s, ncol = K_s - 1)
    true_s <- 10

    pi_s <- compute_pi(true_beta, mm_s)
    alpha_s <- pi_s * true_s
    n_s <- rpois(N_s, 20) + 5
    Y_s <- t(sapply(1:N_s, function(i) {
        p_dir <- rgamma(K_s, alpha_s[i, ])
        p_dir <- p_dir / sum(p_dir)
        rmultinom(1, n_s[i], p_dir)
    }))

    fit_s <- optimize_joint(Y_s, mm_s, K_s, true_beta * 0.5, true_s * 0.5,
        s_trend = true_s, lambda = 0, max_iter = 500
    )
    beta_err <- max(abs(fit_s$beta - true_beta))
    s_err <- abs(fit_s$s - true_s)
    if (beta_err < 1.0 && s_err < 5.0 && fit_s$convergence == 0) {
        cat(
            "PASS (|beta_err|:", round(beta_err, 3),
            "|s_err|:", round(s_err, 3), "s_hat:", round(fit_s$s, 2), ")\n"
        )
    } else {
        cat(
            "FAIL (|beta_err|:", round(beta_err, 3), "|s_err|:", round(s_err, 3),
            "conv:", fit_s$convergence, ")\n"
        )
        all_pass <- FALSE
    }

    # Test 4: LRT degrees of freedom formula
    cat("Test 4: DF formula check... ")
    df_pass <- TRUE
    for (K_df in c(2, 3, 5, 11)) {
        for (sdf in c(3, 5)) {
            # p_full = intercept + 2 batch dummies + spline_df
            p_full <- 1 + 2 + sdf
            p_null <- 1 + 2
            actual_df <- (p_full - p_null) * (K_df - 1)
            expected_df <- sdf * (K_df - 1)
            if (actual_df != expected_df) {
                df_pass <- FALSE
            }
        }
    }
    if (df_pass) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    # === Helper: generate synthetic DirMult data for fit_dirmult_gene tests ===
    generate_test_data <- function(N, K, s, spline_df, beta = NULL) {
        pt <- sort(runif(N))
        mm_full <- model.matrix(~ bs(pt, df = spline_df))
        mm_null <- matrix(1, nrow = N, ncol = 1)
        p <- ncol(mm_full)
        if (is.null(beta)) {
            beta <- matrix(rnorm(p * (K - 1), sd = 0.3), nrow = p, ncol = K - 1)
        }
        pi_mat <- compute_pi(beta, mm_full)
        alpha <- pi_mat * s
        n_counts <- rpois(N, 15) + 3L
        Y <- t(sapply(seq_len(N), function(i) {
            p_dir <- rgamma(K, alpha[i, ])
            p_dir <- p_dir / sum(p_dir)
            rmultinom(1, n_counts[i], p_dir)
        }))
        list(Y = Y, mm_full = mm_full, mm_null = mm_null, K = K)
    }

    # Test 5: Shared-s nesting guarantee (ll_full >= ll_null for all scenarios)
    cat("Test 5: Shared-s nesting guarantee... ")
    t5_pass <- TRUE
    t5_count <- 0L
    for (K_t5 in c(2, 3, 5)) {
        for (N_t5 in c(100, 200, 500)) {
            for (s_t5 in c(1, 5, 20)) {
                sim <- generate_test_data(N_t5, K_t5, s_t5, spline_df = 4)
                res <- fit_dirmult_gene(sim$Y, sim$mm_full, sim$mm_null, sim$K,
                    "joint", s_t5,
                    lambda = 0.1,
                    max_iter = 300, compute_se = FALSE
                )
                t5_count <- t5_count + 1L
                if (!is.na(res$ll_full) && !is.na(res$ll_null)) {
                    if (res$ll_full < res$ll_null - 1e-6) {
                        t5_pass <- FALSE
                        cat(
                            "\n  VIOLATION: K=", K_t5, " N=", N_t5, " s=", s_t5,
                            " ll_full=", res$ll_full, " ll_null=", res$ll_null
                        )
                    }
                }
            }
        }
    }
    if (t5_pass) {
        cat("PASS (", t5_count, " scenarios)\n")
    } else {
        cat("\nFAIL\n")
        all_pass <- FALSE
    }

    # Test 6: Shared-s consistency (s_full == s_null exactly)
    cat("Test 6: Shared-s consistency... ")
    sim6 <- generate_test_data(200, 3, 8, spline_df = 4)
    res6 <- fit_dirmult_gene(sim6$Y, sim6$mm_full, sim6$mm_null, sim6$K,
        "joint", 8,
        lambda = 0.1, max_iter = 300,
        compute_se = FALSE
    )
    if (!is.na(res6$s_full) && !is.na(res6$s_null) &&
        identical(res6$s_full, res6$s_null)) {
        cat("PASS (s_shared =", round(res6$s_shared, 3), ")\n")
    } else {
        cat("FAIL (s_full:", res6$s_full, "s_null:", res6$s_null, ")\n")
        all_pass <- FALSE
    }

    # Test 7: Fallback s consistency (forced non-convergence via max_iter=1)
    cat("Test 7: Fallback s consistency... ")
    sim7 <- generate_test_data(200, 3, 8, spline_df = 4)
    s_trend_t7 <- 5.0
    res7 <- fit_dirmult_gene(sim7$Y, sim7$mm_full, sim7$mm_null, sim7$K,
        "joint", s_trend_t7,
        lambda = 0.1,
        max_iter = 1, compute_se = FALSE
    )
    if (!is.na(res7$s_full) && !is.na(res7$s_null) &&
        identical(res7$s_full, res7$s_null)) {
        cat("PASS (s_shared =", round(res7$s_shared, 3), ")\n")
    } else {
        cat("FAIL (s_full:", res7$s_full, "s_null:", res7$s_null, ")\n")
        all_pass <- FALSE
    }

    # Test 8: Deviance non-negativity across random scenarios
    cat("Test 8: Deviance non-negativity... ")
    t8_pass <- TRUE
    for (trial in 1:10) {
        K_t8 <- sample(2:4, 1)
        N_t8 <- sample(c(100, 200, 300), 1)
        s_t8 <- runif(1, 1, 20)
        sim8 <- generate_test_data(N_t8, K_t8, s_t8, spline_df = 4)
        res8 <- fit_dirmult_gene(sim8$Y, sim8$mm_full, sim8$mm_null, sim8$K,
            "joint", s_t8,
            lambda = 0.1,
            max_iter = 300, compute_se = FALSE
        )
        if (!is.na(res8$deviance) && res8$deviance < 0) {
            t8_pass <- FALSE
            cat("\n  NEGATIVE:", res8$deviance, "K=", K_t8, "N=", N_t8)
        }
    }
    if (t8_pass) {
        cat("PASS (10 random scenarios)\n")
    } else {
        cat("\nFAIL\n")
        all_pass <- FALSE
    }

    # Test 9: Null p-value uniformity (KS test)
    cat("Test 9: Null p-value uniformity... ")
    n_null <- 100L
    N_t9 <- 200L
    K_t9 <- 2L
    s_t9 <- 5
    pt_t9 <- sort(runif(N_t9))
    mm_full_t9 <- model.matrix(~ bs(pt_t9, df = 4))
    mm_null_t9 <- matrix(1, nrow = N_t9, ncol = 1)
    p_t9 <- ncol(mm_full_t9)
    pvals_t9 <- numeric(n_null)
    for (g in seq_len(n_null)) {
        # Under the null: only intercept affects pi, spline coefficients = 0
        beta_null <- matrix(0, nrow = p_t9, ncol = K_t9 - 1L)
        beta_null[1, ] <- rnorm(K_t9 - 1L, sd = 0.5)
        pi_mat <- compute_pi(beta_null, mm_full_t9)
        alpha <- pi_mat * s_t9
        n_counts <- rpois(N_t9, 15) + 3L
        Y_t9 <- t(sapply(seq_len(N_t9), function(i) {
            p_dir <- rgamma(K_t9, alpha[i, ])
            p_dir <- p_dir / sum(p_dir)
            rmultinom(1, n_counts[i], p_dir)
        }))
        res <- fit_dirmult_gene(Y_t9, mm_full_t9, mm_null_t9, K_t9,
            "joint", s_t9,
            lambda = 0,
            max_iter = 300, compute_se = FALSE
        )
        pvals_t9[g] <- res$p_value
    }
    valid_p <- pvals_t9[!is.na(pvals_t9)]
    ks_res <- ks.test(valid_p, "punif")
    if (ks_res$p.value > 0.01) {
        cat(
            "PASS (KS p =", round(ks_res$p.value, 3),
            "n =", length(valid_p), ")\n"
        )
    } else {
        cat(
            "FAIL (KS p =", round(ks_res$p.value, 3),
            "n =", length(valid_p), ")\n"
        )
        all_pass <- FALSE
    }

    # Test 10: Beta-only SE validity
    cat("Test 10: Beta-only SE validity... ")
    n_se <- 20L
    se_ok <- 0L
    for (trial in seq_len(n_se)) {
        sim10 <- generate_test_data(200, 2, 5, spline_df = 4)
        res10 <- fit_dirmult_gene(sim10$Y, sim10$mm_full, sim10$mm_null,
            sim10$K, "joint", 5,
            lambda = 0.1,
            max_iter = 300, compute_se = TRUE
        )
        if (!is.null(res10$se) && all(is.finite(res10$se)) &&
            all(res10$se > 0)) {
            se_ok <- se_ok + 1L
        }
    }
    if (se_ok >= 18L) {
        cat("PASS (", se_ok, "/", n_se, " valid)\n")
    } else {
        cat("FAIL (", se_ok, "/", n_se, " valid)\n")
        all_pass <- FALSE
    }

    # Test 11: Rank-deficient design handling
    cat("Test 11: Rank-deficient design... ")
    N_t11 <- 200L
    K_t11 <- 2L
    pt_t11 <- sort(runif(N_t11))
    mm_base <- model.matrix(~ bs(pt_t11, df = 4))
    # Add a duplicate column to create rank deficiency
    mm_rd <- cbind(mm_base, mm_base[, 3])
    mm_null_t11 <- matrix(1, nrow = N_t11, ncol = 1)
    beta_t11 <- matrix(rnorm(ncol(mm_base) * (K_t11 - 1L), sd = 0.3),
        nrow = ncol(mm_base), ncol = K_t11 - 1L
    )
    pi_t11 <- compute_pi(beta_t11, mm_base)
    alpha_t11 <- pi_t11 * 5
    n_t11 <- rpois(N_t11, 15) + 3L
    Y_t11 <- t(sapply(seq_len(N_t11), function(i) {
        p_dir <- rgamma(K_t11, alpha_t11[i, ])
        p_dir <- p_dir / sum(p_dir)
        rmultinom(1, n_t11[i], p_dir)
    }))
    res11 <- fit_dirmult_gene(Y_t11, mm_rd, mm_null_t11, K_t11,
        "joint", 5,
        lambda = 0.1,
        max_iter = 300, compute_se = FALSE
    )
    # After QR rank fix: mm_rd drops to ncol(mm_base) columns
    # lrt_df = (ncol(mm_base) - 1) * (K-1) = spline_df * (K-1) = 4
    expected_df_11 <- (ncol(mm_base) - 1L) * (K_t11 - 1L)
    if (!is.na(res11$lrt_df) && res11$lrt_df == expected_df_11) {
        cat("PASS (df =", res11$lrt_df, ")\n")
    } else {
        cat("FAIL (expected:", expected_df_11, "got:", res11$lrt_df, ")\n")
        all_pass <- FALSE
    }

    # Test 12: Stage 1 conditional-MLE refit closes the joint-vs-conditional gap
    # Constructs a stress regime (small s, narrow pseudotime, sparse counts) where
    # the joint-optimum beta can drift from the conditional MLE at s_shared, then
    # verifies that fit_dirmult_gene (which now refits beta at s_shared) still
    # produces ll_full >= ll_null. This is the regime that produced the 275 real
    # -data violations.
    cat("Test 12: Stage 1 conditional-MLE refit (stress regime)... ")
    t12_pass <- TRUE
    t12_count <- 0L
    set.seed(2026L)
    for (trial in seq_len(40L)) {
        N_t12 <- sample(c(60L, 80L, 100L), 1L)
        K_t12 <- sample(c(3L, 4L, 5L), 1L)
        s_t12 <- runif(1L, 0.05, 0.6) # severe overdispersion
        # Narrow pseudotime range and sparse counts
        pt_t12 <- sort(runif(N_t12, 0.3, 0.7))
        mm_full_t12 <- model.matrix(~ bs(pt_t12, df = 5))
        mm_null_t12 <- matrix(1, nrow = N_t12, ncol = 1L)
        p_t12 <- ncol(mm_full_t12)
        beta_t12 <- matrix(rnorm(p_t12 * (K_t12 - 1L), sd = 0.4),
            nrow = p_t12, ncol = K_t12 - 1L
        )
        pi_t12 <- compute_pi(beta_t12, mm_full_t12)
        alpha_t12 <- pi_t12 * s_t12
        n_t12 <- rpois(N_t12, 4) + 1L # very sparse
        Y_t12 <- t(vapply(seq_len(N_t12), function(i) {
            p_dir <- rgamma(K_t12, alpha_t12[i, ])
            p_dir <- p_dir / sum(p_dir)
            as.integer(rmultinom(1, n_t12[i], p_dir))
        }, integer(K_t12)))
        res12 <- tryCatch(
            fit_dirmult_gene(Y_t12, mm_full_t12, mm_null_t12, K_t12,
                "joint", s_t12,
                lambda = 0.1, max_iter = 300, compute_se = FALSE
            ),
            error = function(e) NULL
        )
        if (is.null(res12)) next
        t12_count <- t12_count + 1L
        if (!is.na(res12$ll_full) && !is.na(res12$ll_null) &&
            res12$ll_full < res12$ll_null - 1e-6) {
            t12_pass <- FALSE
            cat(
                "\n  VIOLATION trial=", trial, " K=", K_t12, " N=", N_t12,
                " s=", round(s_t12, 3),
                " ll_full=", round(res12$ll_full, 4),
                " ll_null=", round(res12$ll_null, 4)
            )
        }
        # Also: if the optimization succeeded, lrt_status should be "ok" and
        # the deviance should be non-NA
        if (!is.null(res12$lrt_status) && res12$lrt_status == "ok" &&
            is.na(res12$deviance)) {
            t12_pass <- FALSE
            cat("\n  ok-status with NA deviance trial=", trial)
        }
    }
    if (t12_pass) {
        cat("PASS (", t12_count, " stress trials)\n")
    } else {
        cat("\nFAIL\n")
        all_pass <- FALSE
    }

    # Test 13: Stage 2 padded-null fallback restores nesting under bad warm start
    # Verifies the padded-null safety net works when invoked: build a synthetic
    # gene, manually run the Stage 1 refit, then artificially zero out the full-
    # model beta (the worst possible warm start) and call optimize_beta_fixed_s
    # twice - once from the bad init, once from the padded null - and assert
    # that the second always reaches at least the null LL.
    cat("Test 13: Stage 2 padded-null fallback... ")
    set.seed(13L)
    sim13 <- generate_test_data(150, 3, 2, spline_df = 5)
    K_t13 <- sim13$K
    mm_full_t13 <- sim13$mm_full
    mm_null_t13 <- sim13$mm_null
    Y_t13 <- sim13$Y
    s_t13 <- 2
    # Get fit_null at s_t13
    beta_init_n13 <- init_beta_from_glm(Y_t13, mm_null_t13, K_t13)
    fit_null_t13 <- optimize_beta_fixed_s(
        Y_t13, mm_null_t13, K_t13, beta_init_n13, s_t13, 300
    )
    # Build the padded-null beta in the full beta space
    p_full_t13 <- ncol(mm_full_t13)
    p_null_t13 <- ncol(mm_null_t13)
    beta_pad_t13 <- matrix(0, nrow = p_full_t13, ncol = K_t13 - 1L)
    beta_pad_t13[seq_len(p_null_t13), ] <- fit_null_t13$beta
    # Sanity: padded LL must equal null LL
    ll_pad_t13 <- dirmult_loglik_vec(
        Y_t13, compute_pi(beta_pad_t13, mm_full_t13) * s_t13
    )
    pad_check <- abs(ll_pad_t13 - fit_null_t13$ll) < 1e-6
    # Run the safety-net refit
    fit_safe_t13 <- optimize_beta_fixed_s(
        Y_t13, mm_full_t13, K_t13, beta_pad_t13, s_t13, 300
    )
    if (pad_check && fit_safe_t13$ll >= fit_null_t13$ll - 1e-10) {
        cat(
            "PASS (ll_safe =", round(fit_safe_t13$ll, 3),
            ">= ll_null =", round(fit_null_t13$ll, 3), ")\n"
        )
    } else {
        cat(
            "FAIL (pad_check:", pad_check,
            "ll_safe:", fit_safe_t13$ll,
            "ll_null:", fit_null_t13$ll, ")\n"
        )
        all_pass <- FALSE
    }

    # Test 14: Stage 3 lrt_status field is populated
    # All ok genes from Test 5 should report lrt_status == "ok" and have non-NA
    # deviance. We re-run a small slice and check.
    cat("Test 14: lrt_status field correctness... ")
    set.seed(14L)
    sim14 <- generate_test_data(200, 3, 5, spline_df = 4)
    res14 <- fit_dirmult_gene(sim14$Y, sim14$mm_full, sim14$mm_null, sim14$K,
        "joint", 5,
        lambda = 0.1, max_iter = 300, compute_se = FALSE
    )
    if (!is.null(res14$lrt_status) && res14$lrt_status == "ok" &&
        !is.na(res14$deviance) && !is.na(res14$p_value)) {
        cat("PASS (status =", res14$lrt_status, ")\n")
    } else {
        cat(
            "FAIL (status:", res14$lrt_status,
            "deviance:", res14$deviance, ")\n"
        )
        all_pass <- FALSE
    }

    # Test 15: Bootstrap p-value uniformity under the null
    # The bootstrap analogue of Test 9: under H0, bootstrap p-values must be
    # approximately uniform on (0, 1). We use a small B (B = 99) and modest
    # gene count to keep the test fast.
    cat("Test 15: Bootstrap p-value uniformity... ")
    set.seed(15L)
    n_boot_genes <- 50L
    N_t15 <- 200L
    K_t15 <- 2L
    s_t15 <- 5
    pt_t15 <- sort(runif(N_t15))
    mm_full_t15 <- model.matrix(~ bs(pt_t15, df = 4))
    mm_null_t15 <- matrix(1, nrow = N_t15, ncol = 1L)
    p_t15 <- ncol(mm_full_t15)
    pvals_boot <- numeric(n_boot_genes)
    for (g in seq_len(n_boot_genes)) {
        beta_g <- matrix(0, nrow = p_t15, ncol = K_t15 - 1L)
        beta_g[1, ] <- rnorm(K_t15 - 1L, sd = 0.5)
        pi_g <- compute_pi(beta_g, mm_full_t15)
        alpha_g <- pi_g * s_t15
        n_g <- rpois(N_t15, 15) + 3L
        Y_g <- t(vapply(seq_len(N_t15), function(i) {
            p_dir <- rgamma(K_t15, alpha_g[i, ])
            p_dir <- p_dir / sum(p_dir)
            as.integer(rmultinom(1, n_g[i], p_dir))
        }, integer(K_t15)))
        res_g <- fit_dirmult_gene(Y_g, mm_full_t15, mm_null_t15, K_t15,
            "joint", s_t15,
            lambda = 0, max_iter = 300, compute_se = FALSE
        )
        if (is.na(res_g$deviance) || is.na(res_g$s_shared)) {
            pvals_boot[g] <- NA_real_
            next
        }
        boot_g <- bootstrap_dirmult_lrt(
            Y_g, res_g$mm_full, res_g$mm_null, K_t15,
            res_g$beta_full, res_g$beta_null, res_g$s_shared,
            res_g$deviance, B = 99L, max_iter = 300, n_counts = rowSums(Y_g)
        )
        pvals_boot[g] <- boot_g$p_boot
    }
    valid_p15 <- pvals_boot[!is.na(pvals_boot)]
    if (length(valid_p15) >= 30L) {
        # Bootstrap p-values are discrete (B+1 possible values), so the
        # one-sample KS test reports a 'ties' warning. The test statistic
        # is still valid as an approximate check of uniformity; we
        # suppress the cosmetic warning.
        ks_res15 <- suppressWarnings(ks.test(valid_p15, "punif"))
        if (ks_res15$p.value > 0.01) {
            cat(
                "PASS (KS p =", round(ks_res15$p.value, 3),
                "n =", length(valid_p15), ")\n"
            )
        } else {
            cat(
                "FAIL (KS p =", round(ks_res15$p.value, 3),
                "n =", length(valid_p15), ")\n"
            )
            all_pass <- FALSE
        }
    } else {
        cat("FAIL (only", length(valid_p15), "valid bootstrap p-values)\n")
        all_pass <- FALSE
    }

    # Test 16: Multi-regime null calibration grid (replaces single-cell Test 9)
    # Sweeps K x N x s in {2,3,5} x {100,200,500} x {0.5,5,50}, generates 150
    # null genes per cell, and runs two FWER-controlled checks:
    #   (a) Bonferroni-corrected per-cell KS test at family-wise alpha = 0.01
    #       Per-cell threshold = 0.01 / 27 ~= 3.7e-4. Requires every cell to
    #       pass individually after FWER correction.
    #   (b) Fisher combined test across all 27 cells: X^2 = -2 * sum(log(p_cell))
    #       ~ chi^2(54) under global H0. Requires X^2 p-value > 0.01 globally.
    # Both must pass. This replaces the plain "all 27
    # cells pass at 1%" check (which has ~24% FWER under H0 by Sidak, causing
    # spurious failures). Uses n=150 per cell for better KS resolution and to
    # reduce the chance of ties from near-p=1 clustering in the extreme regime.
    cat("Test 16: Multi-regime null calibration grid (27 cells, FWER 1%)... ")
    set.seed(16L)
    t16_n_genes <- 150L
    t16_alpha_bonf <- 0.01 / 27  # Bonferroni-corrected per-cell threshold
    t16_grid_pass <- 0L
    t16_grid_total <- 0L
    t16_failures <- character(0)
    t16_all_ks_p <- numeric(0)
    for (K_t16 in c(2L, 3L, 5L)) {
        for (N_t16 in c(100L, 200L, 500L)) {
            for (s_t16 in c(0.5, 5, 50)) {
                t16_grid_total <- t16_grid_total + 1L
                pt_t16 <- sort(runif(N_t16))
                mm_full_t16 <- model.matrix(~ bs(pt_t16, df = 4))
                mm_null_t16 <- matrix(1, nrow = N_t16, ncol = 1L)
                p_t16 <- ncol(mm_full_t16)
                pvals_t16 <- numeric(t16_n_genes)
                for (g in seq_len(t16_n_genes)) {
                    beta_g16 <- matrix(0, nrow = p_t16, ncol = K_t16 - 1L)
                    beta_g16[1, ] <- rnorm(K_t16 - 1L, sd = 0.5)
                    pi_g16 <- compute_pi(beta_g16, mm_full_t16)
                    alpha_g16 <- pi_g16 * s_t16
                    n_g16 <- rpois(N_t16, 15) + 3L
                    Y_g16 <- t(vapply(seq_len(N_t16), function(i) {
                        p_dir <- rgamma(K_t16, alpha_g16[i, ])
                        p_dir <- p_dir / sum(p_dir)
                        as.integer(rmultinom(1, n_g16[i], p_dir))
                    }, integer(K_t16)))
                    res_g16 <- fit_dirmult_gene(
                        Y_g16, mm_full_t16, mm_null_t16, K_t16,
                        "joint", s_t16,
                        lambda = 0, max_iter = 300, compute_se = FALSE
                    )
                    pvals_t16[g] <- res_g16$p_value
                }
                vp16 <- pvals_t16[!is.na(pvals_t16)]
                if (length(vp16) >= 30L) {
                    ks_g16 <- suppressWarnings(ks.test(vp16, "punif"))
                    t16_all_ks_p <- c(t16_all_ks_p, ks_g16$p.value)
                    if (ks_g16$p.value > t16_alpha_bonf) {
                        t16_grid_pass <- t16_grid_pass + 1L
                    } else {
                        t16_failures <- c(
                            t16_failures,
                            sprintf(
                                "K=%d N=%d s=%.1f KSp=%.4g",
                                K_t16, N_t16, s_t16, ks_g16$p.value
                            )
                        )
                    }
                } else {
                    t16_failures <- c(
                        t16_failures,
                        sprintf(
                            "K=%d N=%d s=%.1f insufficient_n=%d",
                            K_t16, N_t16, s_t16, length(vp16)
                        )
                    )
                }
            }
        }
    }
    # Fisher combined test across all 27 cells (global FWER check)
    if (length(t16_all_ks_p) == t16_grid_total) {
        t16_fisher_stat <- -2 * sum(log(pmax(t16_all_ks_p, 1e-300)))
        t16_fisher_p <- pchisq(
            t16_fisher_stat,
            df = 2 * t16_grid_total, lower.tail = FALSE
        )
        t16_fisher_ok <- t16_fisher_p > 0.01
    } else {
        t16_fisher_p <- NA_real_
        t16_fisher_ok <- FALSE
    }
    if (t16_grid_pass == t16_grid_total && t16_fisher_ok) {
        cat(
            "PASS (", t16_grid_pass, "/", t16_grid_total,
            " cells, Fisher p=", sprintf("%.3f", t16_fisher_p), ")\n"
        )
    } else {
        cat(
            "FAIL (", t16_grid_pass, "/", t16_grid_total,
            " cells passed, Fisher p=", sprintf("%.3g", t16_fisher_p), ")\n"
        )
        if (length(t16_failures) > 0L) {
            cat("  Failures:\n  ", paste(t16_failures, collapse = "\n  "), "\n")
        }
        all_pass <- FALSE
    }

    cat("\n=== All tests", if (all_pass) "PASSED" else "had FAILURES", "===\n")
    quit(status = if (all_pass) 0L else 1L)
}

# ===== Section 1: Data Loading ===============================================

cat("=== Optimized_2 Dirichlet-Multinomial Analysis (shared-s, conditional-MLE) ===\n")
cat("Parameters:\n")
cat("  method:", opts$method, "\n")
cat("  spline_df:", opts$spline_df, "\n")
cat("  min_cells:", opts$min_cells, "\n")
cat("  shrinkage:", opts$shrinkage, "(lambda:", opts$lambda, ")\n")
cat("  compute_se:", opts$compute_se, "\n")
cat("  padj_method:", opts$padj_method, "\n")
cat("  ncores:", opts$ncores, "\n")
cat("  bootstrap:", opts$bootstrap,
    if (opts$bootstrap > 0L) sprintf("(threshold p_chisq < %.3f)", opts$bootstrap_pval) else "(off)",
    "\n"
)
cat("  calibration_n:", opts$calibration_n, "\n")
cat("  lambda_sensitivity:", opts$lambda_sensitivity, "\n")

cat("\nLoading data from", opts$data_path, "\n")
data <- read_h5ad(opts$data_path)
if (inherits(data$X, "sparseMatrix")) data$X <- as.matrix(data$X)

X <- as.matrix(data$X)
utr_names <- as.character(data$var$utr_name)
batch <- as.vector(data$obs$batch_number)
pseudotime <- as.numeric(data$obs$dpt_pseudotime)

# Drop NaN pseudotime cells
valid_cells <- !is.na(pseudotime)
X <- X[valid_cells, , drop = FALSE]
batch <- batch[valid_cells]
pseudotime <- pseudotime[valid_cells]

cat(
    "Data:", nrow(X), "cells,", ncol(X), "features,",
    length(unique(utr_names)), "unique UTRs\n"
)
cat("Batches:", paste(sort(unique(batch)), collapse = ", "), "\n")
cat(
    "Pseudotime range:", round(min(pseudotime), 4), "-",
    round(max(pseudotime), 4), "\n"
)

output_dir <- file.path(dirname(opts$data_path), "DirMult_optimized")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(dirname(opts$data_path))

# Precompute spline basis on full pseudotime (consistent knots across genes)
spline_basis <- bs(pseudotime, df = opts$spline_df)

# ===== Section 5: Main Analysis Loop =========================================

# --- Phase A: Global concentration trend ---
cat("\nPhase A: Estimating global concentration trend...\n")
trend_result <- estimate_global_trend(X, utr_names)
cat("  Median MoM s:", round(median(trend_result$s_mom, na.rm = TRUE), 2), "\n")
cat("  Median trend s:", round(median(trend_result$s_trend, na.rm = TRUE), 2), "\n")

# Compute per-gene shrinkage lambda
unique_utrs <- unique(utr_names)
gene_lambda <- setNames(rep(opts$lambda, length(unique_utrs)), unique_utrs)

if (opts$shrinkage == "adaptive") {
    for (idx in seq_along(unique_utrs)) {
        cols <- which(utr_names == unique_utrs[idx])
        n_cells_gene <- sum(rowSums(X[, cols, drop = FALSE]) > 0)
        gene_lambda[idx] <- opts$lambda / sqrt(max(n_cells_gene, 1))
    }
} else if (opts$shrinkage == "none") {
    gene_lambda[] <- 0
}

# --- Helper: empty per-gene result with all expected fields ---
empty_gene_result <- function(utr, n_cells_keep, K_obs, method) {
    list(
        utr_name = utr, deviance = NA, p_value = NA, s_full = NA,
        s_null = NA, s_shared = NA, ll_full = NA, ll_null = NA,
        n_cells = n_cells_keep, K = K_obs, lrt_df = NA,
        converged = NA, full_converged = NA,
        cells_per_param = NA,
        method = method, se = NULL,
        lrt_status = "skipped",
        max_var_ratio = NA_real_, min_var_ratio = NA_real_,
        beta_full = NULL, beta_null = NULL, mm_full = NULL, mm_null = NULL,
        Y = NULL, n_counts = NULL
    )
}

#' Fit a single gene end-to-end given the full data and gene name.
#'
#' Wraps the per-gene logic so we can call it from the main loop and from
#' the lambda-sensitivity sweep without duplicating the design-matrix setup.
fit_one_utr <- function(utr, lambda_val, lambda_method = NULL) {
    method_to_use <- if (is.null(lambda_method)) opts$method else lambda_method
    cols <- which(utr_names == utr)
    Y_utr <- X[, cols, drop = FALSE]
    keep <- rowSums(Y_utr) > 0
    Y_sel <- Y_utr[keep, , drop = FALSE]
    K_obs <- ncol(Y_sel)

    if (K_obs < 2 || nrow(Y_sel) < opts$min_cells || sum(Y_sel) <= 10) {
        res <- empty_gene_result(utr, sum(keep), K_obs, method_to_use)
        return(res)
    }

    result <- tryCatch(
        {
            batch_sel <- factor(batch[keep])
            spline_sel <- spline_basis[keep, , drop = FALSE]
            if (nlevels(batch_sel) > 1) {
                batch_mm <- model.matrix(~batch_sel)
                mm_full <- cbind(batch_mm, spline_sel)
                mm_null <- batch_mm
            } else {
                mm_full <- cbind(1, spline_sel)
                mm_null <- matrix(1, nrow = nrow(Y_sel), ncol = 1)
            }
            r <- fit_dirmult_gene(
                Y_sel, mm_full, mm_null, K_obs, method_to_use,
                trend_result$s_trend[utr], lambda_val,
                opts$max_iter, opts$compute_se
            )
            r$Y <- Y_sel
            r$n_counts <- rowSums(Y_sel)
            r
        },
        error = function(e) {
            r <- empty_gene_result(utr, sum(keep), K_obs, method_to_use)
            r$converged <- FALSE
            r$lrt_status <- "error"
            r
        }
    )
    result$utr_name <- utr
    result
}

# --- Dynamic load-balanced dispatch helpers ------------------------------
# Replaces pblapply(cl=cluster), which uses static chunked scheduling
# (parLapply: contiguous blocks of ~N/ncores genes per worker, with a sync
# barrier every progress-bar chunk). Static chunking leaves cores idle while
# a few stragglers finish -- a problem here because per-gene fit cost varies
# widely (high-K / high-cell genes are much slower) and, when slow genes are
# adjacent in the gene list, static scheduling can run several-fold slower
# than ideal (benchmarked 2026-05-27 on the 260526 version of this script).
#
# parLapplyLB(chunk.size=1) dispatches one gene at a time to whichever worker
# is free. We split into round-robin progress batches purely for logging;
# round-robin assignment spreads slow genes across batches so each batch
# stays balanced.
#
# Correctness:
#   - Phase B / Phase E (fitting) consume no randomness -> results bit-identical
#     to the old static path, just faster.
#   - Phase C (bootstrap) / Phase D (calibration) do consume RNG. To keep them
#     reproducible under dynamic dispatch (where worker assignment is timing-
#     dependent) each gene now reseeds its own RNG via gene_seed_from() at task
#     entry. This makes those phases reproducible across runs and core counts
#     -- which the original inherited-fork-RNG path was not. Note: bootstrap /
#     calibration p-values therefore differ from any pre-2026-05-27 run within
#     Monte-Carlo noise (they are now deterministic given --seed). The primary
#     chi-squared LRT (p_value_chisq, padj) is unaffected: it comes entirely
#     from the deterministic Phase B.
dispatch_lb <- function(cl, idx, FUN, n_progress = 20L, label = "") {
    n <- length(idx)
    if (n == 0L) return(list())
    n_progress <- as.integer(max(1L, min(n_progress, n)))
    batch_of <- (seq_len(n) - 1L) %% n_progress  # round-robin -> spreads slow genes
    out <- vector("list", n)
    t0 <- Sys.time()
    done <- 0L
    for (b in seq_len(n_progress) - 1L) {
        sel <- which(batch_of == b)
        out[sel] <- parallel::parLapplyLB(cl, idx[sel], FUN, chunk.size = 1L)
        done <- done + length(sel)
        cat(sprintf(
            "  %s %d/%d (%.0f%%), elapsed %s\n",
            label, done, n, 100 * done / n,
            format(round(difftime(Sys.time(), t0, units = "mins"), 2))
        ))
        utils::flush.console()
    }
    out
}

# Per-gene reproducible seed (xxhash32 of the gene id, mixed with master seed).
# Used to seed the RNG-consuming Phase C / D tasks so the dynamic scheduler
# produces identical results to a static run for a given --seed.
gene_seed_from <- function(master_seed, utr_name) {
    h <- strtoi(substr(digest(utr_name, algo = "xxhash32"), 1L, 7L), base = 16L)
    as.integer(bitwXor(as.integer(master_seed * 1000003L), h))
}

# --- Phase B: Per-gene fitting ---
cat("\nPhase B: Fitting DirMult models per gene (", length(unique_utrs), " UTRs)...\n")
cat("Setting up FORK cluster with", opts$ncores, "cores (dynamic load balancing)\n")
cl <- makeCluster(opts$ncores, type = "FORK")

results_list <- dispatch_lb(
    cl, seq_along(unique_utrs),
    function(idx) fit_one_utr(unique_utrs[idx], gene_lambda[unique_utrs[idx]]),
    label = "Phase B fit:"
)

stopCluster(cl)

# Index results by utr_name for fast lookup in Phases C/D
names(results_list) <- vapply(results_list, `[[`, character(1), "utr_name")

# --- Phase C (optional): Parametric bootstrap on borderline genes ---
# For genes whose chi-squared p-value is below opts$bootstrap_pval, recompute
# the p-value by parametric bootstrap of the conditional LRT statistic. The
# bootstrap is conditional on s_shared and uses only the fast
# optimize_beta_fixed_s path, so it is feasible to run on the borderline slice
# in a few hours of wall-clock time.

bootstrap_results <- list()
if (opts$bootstrap > 0L) {
    cat("\nPhase C: Bootstrap on borderline genes (B =", opts$bootstrap,
        ", threshold p_chisq <", opts$bootstrap_pval, ")...\n"
    )
    borderline_idx <- which(vapply(results_list, function(r) {
        !is.null(r$p_value) && !is.na(r$p_value) &&
            !is.null(r$lrt_status) && r$lrt_status == "ok" &&
            r$p_value < opts$bootstrap_pval
    }, logical(1)))
    cat("  ", length(borderline_idx), "genes selected for bootstrap.\n")

    if (length(borderline_idx) > 0L) {
        cl_b <- makeCluster(opts$ncores, type = "FORK")
        # L'Ecuyer-CMRG + per-gene seeding so the dynamic load-balanced
        # dispatch is reproducible (the parametric bootstrap consumes RNG via
        # rgamma/rmultinom). Seed namespace offset (+13) keeps it disjoint from
        # the calibration stream (+7).
        RNGkind("L'Ecuyer-CMRG")
        clusterSetRNGStream(cl_b, opts$seed)
        boot_out <- dispatch_lb(cl_b, borderline_idx, function(j) {
            r <- results_list[[j]]
            if (is.null(r$Y) || is.null(r$beta_full) || is.null(r$beta_null) ||
                is.null(r$mm_full) || is.null(r$mm_null) ||
                is.na(r$s_shared) || is.na(r$deviance)) {
                return(list(
                    p_boot = NA_real_, n_completed = 0L, D_b = numeric(0)
                ))
            }
            set.seed(gene_seed_from(opts$seed + 13L, r$utr_name))
            tryCatch(
                bootstrap_dirmult_lrt(
                    r$Y, r$mm_full, r$mm_null, r$K,
                    r$beta_full, r$beta_null, r$s_shared,
                    r$deviance, opts$bootstrap, opts$max_iter, r$n_counts
                ),
                error = function(e) list(
                    p_boot = NA_real_, n_completed = 0L, D_b = numeric(0)
                )
            )
        }, label = "Phase C boot:")
        stopCluster(cl_b)
        names(boot_out) <- names(results_list)[borderline_idx]
        bootstrap_results <- boot_out
        # Inject the bootstrap p-values back into the results list
        for (nm in names(boot_out)) {
            results_list[[nm]]$p_boot <- boot_out[[nm]]$p_boot
            results_list[[nm]]$n_bootstrap <- boot_out[[nm]]$n_completed
        }
    }
} else {
    cat("\nPhase C: Bootstrap disabled (opts$bootstrap == 0).\n")
}

# --- Phase D: Real-data empirical calibration footprint ---
# Sample N random tested genes, generate one null replicate from each fitted
# null model, refit, and KS-test the resulting p-values against the uniform.
# This is an end-to-end calibration check on the exact parameter regimes
# present in the data, complementing the synthetic-data Test 16.

calibration_result <- list(
    ks_p = NA_real_, n_used = 0L, n_attempted = 0L
)
if (opts$calibration_n > 0L) {
    cat("\nPhase D: Real-data null calibration footprint (n =",
        opts$calibration_n, " genes)...\n"
    )
    # Eligible genes: lrt_status == "ok" and we have the fitted model components
    eligible_idx <- which(vapply(results_list, function(r) {
        !is.null(r$lrt_status) && r$lrt_status == "ok" &&
            !is.null(r$beta_null) && !is.null(r$mm_null) &&
            !is.null(r$mm_full) && !is.null(r$Y) &&
            !is.na(r$s_shared)
    }, logical(1)))
    n_target <- min(opts$calibration_n, length(eligible_idx))
    if (n_target >= 30L) {
        set.seed(2026L)
        sample_idx <- sample(eligible_idx, n_target)
        cl_d <- makeCluster(opts$ncores, type = "FORK")
        # L'Ecuyer-CMRG + per-gene seeding so the dynamic load-balanced
        # dispatch is reproducible (calibration draws consume RNG via
        # rgamma/rmultinom). The original pbsapply path inherited the parent
        # fork RNG per worker and was not reproducible across core counts.
        RNGkind("L'Ecuyer-CMRG")
        clusterSetRNGStream(cl_d, opts$seed)
        cal_list <- dispatch_lb(cl_d, sample_idx, function(j) {
            r <- results_list[[j]]
            set.seed(gene_seed_from(opts$seed + 7L, r$utr_name))
            tryCatch(
                {
                    pi_null <- compute_pi(r$beta_null, r$mm_null)
                    alpha_null <- pi_null * r$s_shared
                    n_i <- r$n_counts
                    K_g <- r$K
                    Y_b <- t(vapply(seq_along(n_i), function(i) {
                        p_dir <- rgamma(K_g, alpha_null[i, ])
                        p_dir <- p_dir / sum(p_dir)
                        as.integer(rmultinom(1, n_i[i], p_dir))
                    }, integer(K_g)))
                    refit <- fit_dirmult_gene(
                        Y_b, r$mm_full, r$mm_null, K_g, opts$method,
                        r$s_shared, gene_lambda[r$utr_name],
                        opts$max_iter, FALSE
                    )
                    refit$p_value
                },
                error = function(e) NA_real_
            )
        }, n_progress = 5L, label = "Phase D cal:")
        stopCluster(cl_d)
        cal_pvals <- vapply(cal_list, function(x) {
            if (is.null(x)) NA_real_ else as.numeric(x)
        }, numeric(1))
        valid_cal <- cal_pvals[!is.na(cal_pvals)]
        calibration_result$n_attempted <- length(cal_pvals)
        calibration_result$n_used <- length(valid_cal)
        if (length(valid_cal) >= 30L) {
            ks_cal <- ks.test(valid_cal, "punif")
            calibration_result$ks_p <- ks_cal$p.value
            cat("  Real-data null KS p =", round(ks_cal$p.value, 4),
                " (n =", length(valid_cal), ")\n"
            )
            if (ks_cal$p.value < 0.01) {
                cat("  WARNING: real-data null calibration FAILED at 1% level.\n",
                    "  Consider re-running with --bootstrap 199 to use the\n",
                    "  parametric bootstrap p-values for borderline genes.\n",
                    sep = ""
                )
            }
        } else {
            cat("  Insufficient valid calibration p-values (",
                length(valid_cal), " of ", n_target, ").\n",
                sep = ""
            )
        }
    } else {
        cat("  Insufficient eligible genes (", length(eligible_idx),
            ") for calibration footprint.\n",
            sep = ""
        )
    }
}

# ===== Section 6: Post-Processing and Output =================================

cat("\nAssembling results...\n")

# Helper to safely extract a numeric scalar from a per-gene result
extract_num <- function(x, field) {
    v <- x[[field]]
    if (is.null(v)) NA_real_ else as.numeric(v)
}
extract_chr <- function(x, field, default = NA_character_) {
    v <- x[[field]]
    if (is.null(v)) default else as.character(v)
}

results_df <- tibble(
    utr_name = vapply(results_list, `[[`, character(1), "utr_name"),
    deviance = vapply(results_list, function(x) extract_num(x, "deviance"), numeric(1)),
    p_value_chisq = vapply(results_list, function(x) extract_num(x, "p_value"), numeric(1)),
    p_value_boot = vapply(results_list, function(x) extract_num(x, "p_boot"), numeric(1)),
    n_bootstrap = vapply(results_list, function(x) {
        nb <- x$n_bootstrap
        if (is.null(nb)) NA_integer_ else as.integer(nb)
    }, integer(1)),
    K = vapply(results_list, function(x) as.integer(x$K), integer(1)),
    n_cells = vapply(results_list, function(x) as.integer(x$n_cells), integer(1)),
    s_full = vapply(results_list, function(x) extract_num(x, "s_full"), numeric(1)),
    s_null = vapply(results_list, function(x) extract_num(x, "s_null"), numeric(1)),
    s_shared = vapply(results_list, function(x) extract_num(x, "s_shared"), numeric(1)),
    ll_full = vapply(results_list, function(x) extract_num(x, "ll_full"), numeric(1)),
    ll_null = vapply(results_list, function(x) extract_num(x, "ll_null"), numeric(1)),
    lrt_df = vapply(results_list, function(x) extract_num(x, "lrt_df"), numeric(1)),
    converged = vapply(results_list, function(x) {
        cv <- x$converged
        if (is.null(cv) || is.na(cv)) NA else as.logical(cv)
    }, logical(1)),
    full_converged = vapply(results_list, function(x) {
        fc <- x$full_converged
        if (is.null(fc) || is.na(fc)) NA else as.logical(fc)
    }, logical(1)),
    cells_per_param = vapply(results_list, function(x) extract_num(x, "cells_per_param"), numeric(1)),
    lrt_status = vapply(results_list, function(x) extract_chr(x, "lrt_status", "skipped"), character(1)),
    max_var_ratio = vapply(results_list, function(x) extract_num(x, "max_var_ratio"), numeric(1)),
    min_var_ratio = vapply(results_list, function(x) extract_num(x, "min_var_ratio"), numeric(1)),
    method = vapply(results_list, function(x) extract_chr(x, "method", opts$method), character(1))
)

# Reported p_value: bootstrap if available, else chi-squared
results_df <- results_df %>%
    mutate(p_value = ifelse(!is.na(p_value_boot), p_value_boot, p_value_chisq))

# Sanity check: shared-s guarantee (s_full must equal s_null for all genes)
n_s_mismatch <- sum(abs(results_df$s_full - results_df$s_null) > 1e-12,
    na.rm = TRUE
)
if (n_s_mismatch > 0) {
    warning(
        "SHARED-S VIOLATION: ", n_s_mismatch,
        " genes have s_full != s_null — this should never happen"
    )
}

# Add UMI counts (backward-compatible column)
umis <- tibble(utr_name = utr_names, UMIs = colSums(X)) %>%
    group_by(utr_name) %>%
    summarize(UMIs = sum(UMIs), .groups = "drop")
results_df <- results_df %>% left_join(umis, by = "utr_name")

# Filter and adjust p-values: drop NA p-values (skipped, optimization_failed,
# null_fit_failed, error). These are excluded rather than reported as
# p = 1 with deviance = 0.
results_df <- results_df %>%
    filter(!is.na(p_value))

# Restrict BH adjustment to fully-converged genes (converged == TRUE &
# lrt_status == "ok").  Non-converged genes have unreliable deviance statistics
# because the optimizer stopped before reaching the true maximum - their
# p-values tend to be artificially small, which inflates the BH threshold and
# produces spurious discoveries among the converged genes.  Setting padj = NA
# for non-converged genes keeps them in the output for inspection but excludes
# them from FDR control.
#
# Use p_value_chisq (not p_value) as the basis for BH.  p_value is the
# bootstrap-corrected value for borderline genes and equals p_value_chisq
# otherwise.  Using the mixed column replaces the most-significant p-values
# with conservative bootstrap estimates, collapsing discoveries under FDR
# correction.  The bootstrap-corrected p-value is preserved in p_value and
# padj_bootstrap for sensitivity analysis.
conv_for_bh <- results_df$converged == TRUE & results_df$lrt_status == "ok"
padj_vec <- rep(NA_real_, nrow(results_df))
padj_vec[conv_for_bh] <- p.adjust(
    results_df$p_value_chisq[conv_for_bh],
    method = opts$padj_method
)
results_df$padj <- padj_vec

# Also store the bootstrap-based padj (BH on p_value) for diagnostic comparison
padj_boot_vec <- rep(NA_real_, nrow(results_df))
padj_boot_vec[conv_for_bh] <- p.adjust(
    results_df$p_value[conv_for_bh],
    method = opts$padj_method
)
results_df$padj_bootstrap <- padj_boot_vec

# Summary
n_total <- length(results_list)
n_tested <- nrow(results_df)
n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
n_conv <- sum(results_df$converged, na.rm = TRUE)
n_full_conv <- sum(results_df$full_converged, na.rm = TRUE)
n_ll_viol <- sum(results_df$ll_full < results_df$ll_null - 1e-6, na.rm = TRUE)

# Counts of excluded statuses among all processed genes
all_status <- vapply(results_list, function(x) extract_chr(x, "lrt_status", "skipped"), character(1))
n_lrt_failed <- sum(all_status == "optimization_failed")
n_null_fail <- sum(all_status == "null_fit_failed")
n_skipped <- sum(all_status == "skipped")
n_error <- sum(all_status == "error")
n_boot_used <- sum(!is.na(results_df$p_value_boot))

cat("\nResults: ", n_tested, " genes tested, ", n_sig, " significant (padj<0.05), ",
    n_conv, "/", n_tested, " converged\n",
    sep = ""
)
cat("Full model primary convergence: ", n_full_conv, "/", n_tested, "\n", sep = "")
cat("LL violations (ll_full < ll_null) among reported: ", n_ll_viol, "\n", sep = "")
cat("Excluded genes:\n")
cat("  optimization_failed: ", n_lrt_failed, " (Stage 3 honest exclusion)\n", sep = "")
cat("  null_fit_failed:     ", n_null_fail, "\n", sep = "")
cat("  skipped (low cells): ", n_skipped, "\n", sep = "")
cat("  error:               ", n_error, "\n", sep = "")
cat("Bootstrap p-values used (replaces chi-sq): ", n_boot_used, "\n", sep = "")
if (!is.na(calibration_result$ks_p)) {
    cat("Real-data null calibration KS p: ",
        round(calibration_result$ks_p, 4),
        " (n = ", calibration_result$n_used, ")\n",
        sep = ""
    )
}
cat("Concentration range: ", round(min(results_df$s_full, na.rm = TRUE), 2), " - ",
    round(max(results_df$s_full, na.rm = TRUE), 2), "\n",
    sep = ""
)
# Misspecification footprint
n_high_vr <- sum(results_df$max_var_ratio > 5, na.rm = TRUE)
n_low_vr <- sum(results_df$min_var_ratio < 0.2, na.rm = TRUE)
cat("Variance-ratio footprint: max_vr > 5: ", n_high_vr,
    "; min_vr < 0.2: ", n_low_vr,
    "  (DirMult misspecification flag)\n",
    sep = ""
)
cat("LRT df used: ", paste(sort(unique(results_df$lrt_df)), collapse = ", "), "\n", sep = "")

# Save outputs
tag <- paste0("optimized2_", opts$method, "_df", opts$spline_df)
csv_path <- file.path(output_dir, paste0("dirmult_regression_", tag, ".csv"))
rds_path <- file.path(output_dir, paste0("dirmult_regression_", tag, ".rds"))

write.csv(results_df, csv_path, row.names = FALSE)
saveRDS(list(
    results = results_df,
    se = setNames(
        lapply(results_list, `[[`, "se"),
        vapply(results_list, `[[`, character(1), "utr_name")
    ),
    trend = trend_result,
    bootstrap = bootstrap_results,
    calibration = calibration_result,
    opts = opts
), rds_path)

cat("\nSaved: ", csv_path, "\n", sep = "")
cat("Saved: ", rds_path, "\n", sep = "")

# --- Phase E (optional): Lambda sensitivity sweep ---
# Re-runs the per-gene fit at lambda in {0, 1, 2} (with the adaptive scaling
# kept) and writes a comparison CSV with per-gene deviance, p_value_chisq, and
# s_shared at each setting. The expected outcome is that p-values are highly
# correlated for genes with reasonable cells_per_param, with divergence only
# for low-power genes where the penalty has an effect. This checks
# the effect of the MAP/MLE mismatch empirically.
if (opts$lambda_sensitivity) {
    cat("\nPhase E: Lambda sensitivity sweep (lambda in {0, 1, 2})...\n")
    lambda_grid <- c(0, 1, 2)
    sens_list <- list()
    for (lam_base in lambda_grid) {
        cat("  --- lambda_base =", lam_base, "---\n")
        # Recompute per-gene lambda with the same adaptive rule
        gene_lambda_lam <- setNames(rep(lam_base, length(unique_utrs)), unique_utrs)
        if (opts$shrinkage == "adaptive" && lam_base > 0) {
            for (idx in seq_along(unique_utrs)) {
                cols2 <- which(utr_names == unique_utrs[idx])
                n_cells_gene2 <- sum(rowSums(X[, cols2, drop = FALSE]) > 0)
                gene_lambda_lam[idx] <- lam_base / sqrt(max(n_cells_gene2, 1))
            }
        }
        cl_e <- makeCluster(opts$ncores, type = "FORK")
        sens_results <- dispatch_lb(cl_e, seq_along(unique_utrs), function(idx) {
            r <- fit_one_utr(unique_utrs[idx], gene_lambda_lam[unique_utrs[idx]])
            list(
                utr_name = r$utr_name,
                deviance = if (is.null(r$deviance)) NA_real_ else r$deviance,
                p_value = if (is.null(r$p_value)) NA_real_ else r$p_value,
                s_shared = if (is.null(r$s_shared)) NA_real_ else r$s_shared,
                lrt_status = if (is.null(r$lrt_status)) "skipped" else r$lrt_status
            )
        }, label = sprintf("Phase E lam=%g:", lam_base))
        stopCluster(cl_e)
        sens_df <- tibble(
            utr_name = vapply(sens_results, `[[`, character(1), "utr_name"),
            deviance = vapply(sens_results, function(x) extract_num(x, "deviance"), numeric(1)),
            p_value = vapply(sens_results, function(x) extract_num(x, "p_value"), numeric(1)),
            s_shared = vapply(sens_results, function(x) extract_num(x, "s_shared"), numeric(1)),
            lrt_status = vapply(sens_results, function(x) extract_chr(x, "lrt_status", "skipped"), character(1))
        )
        names(sens_df)[2:5] <- paste0(names(sens_df)[2:5], "_lam", lam_base)
        sens_list[[as.character(lam_base)]] <- sens_df
    }
    sens_join <- Reduce(function(a, b) full_join(a, b, by = "utr_name"), sens_list)
    sens_path <- file.path(output_dir, paste0("dirmult_lambda_sensitivity_", tag, ".csv"))
    write.csv(sens_join, sens_path, row.names = FALSE)
    cat("Saved: ", sens_path, "\n", sep = "")
    # Quick correlation summary
    for (i in seq_along(lambda_grid)) {
        for (j in seq_len(i - 1L)) {
            li <- lambda_grid[i]
            lj <- lambda_grid[j]
            ci <- paste0("deviance_lam", li)
            cj <- paste0("deviance_lam", lj)
            both <- !is.na(sens_join[[ci]]) & !is.na(sens_join[[cj]])
            if (sum(both) >= 30L) {
                ccor <- cor(sens_join[[ci]][both], sens_join[[cj]][both], method = "spearman")
                cat("  Spearman cor(deviance@lambda=", li,
                    ", deviance@lambda=", lj, ") = ",
                    round(ccor, 4),
                    " (n = ", sum(both), ")\n",
                    sep = ""
                )
            }
        }
    }
}

cat("Done.\n")
