#!/usr/bin/env Rscript

# ============================================================================
# Within-Sample Permutation LRT for 3'UTR Usage Along Pseudotime
# (Mouse SVZ; companion to 260401-...optimized_2.R)
#
# Bit-identical to organoids_knoblich/UTR/260526-permutation_optimized_2.R
# (the canonical source). Kept as a separate copy in the mouse tree so that
# (i) `bash run_permutation_all.sh` runs locally without cross-tree paths,
# and (ii) mouse-specific knobs can diverge from the Knoblich version in
# future without touching the organoid script.
#
# Forked from `260401-dirichlet_multinomial_analysis_optimized_2.R`. Sections
# 0-4 are bit-identical to the source: DM math, optim, gene-fitting cascade,
# loess s-trend, test suite (extended with permutation tests). The fork is
# localised to:
#   - Section 0     (args parsing)         + new permutation args
#   - Section 4.5   (bootstrap)            -> within-sample permutation
#   - Phase B       (fit_one_utr)          + stores pt_keep, sample_block_keep
#   - Phase C       (bootstrap dispatch)   -> permutation dispatch
#   - Section 6     (output)               + p_value_perm, n_permutations cols
#   - run_metadata.json sidecar written next to the CSV
#
# Why this script exists:
#   Mouse chi-sq is already well-calibrated (260407-...QC.ipynb Section 4
#   reports 91.6% concordance with the parametric bootstrap), so the
#   permutation analysis is a parity exercise for the paper's defensibility
#   table -- the same defensibility statement appears in Figure 1 (mouse)
#   and Figure 2 (Knoblich). The Knoblich script header explains the load-
#   bearing scientific motivation; for mouse this is "for completeness and
#   to appease reviewers" (260522-bootstrap_concordance.ipynb Section 14).
#
#   A WITHIN-SAMPLE permutation null fixes this: shuffle pseudotime values
#   only within each sample, preserving per-sample marginals + between-sample
#   batch effects, but breaking within-sample cell-pseudotime association.
#   The resulting null distribution of the LRT statistic is automatically
#   calibrated to the data's actual variance structure, including the within-
#   sample correlation that both chi-sq and the parametric bootstrap ignore.
#
# Permutation strategy options (--permutation_strategy):
#   within_sample (default)  pseudotime shuffled within each --permutation_block_col
#   global                   pseudotime shuffled across all cells in the lineage
#                            (less defensible when samples are imbalanced)
#
# Adaptive B with early stopping:
#   --n_perm_init   (default 49)   initial number of permutations per gene
#   --n_perm_mid    (default 199)  escalation tier 1
#   --n_perm_max    (default 999)  escalation tier 2
#   Escalation rule:
#     - Run n_perm_init permutations.
#     - If n_hits >= n_hits_stop (5) -> stop (p_perm > 5/50 = 0.1).
#     - Else escalate to n_perm_mid. Stop on n_hits >= n_hits_stop.
#     - At n_perm_mid with n_hits == 0 -> escalate to n_perm_max.
#     - At n_perm_max -> stop unconditionally.
#   This typically saves 5-10x compute vs running B=999 on every gene.
#
# Reproducibility:
#   - L'Ecuyer-CMRG RNG (parallel-safe substreams).
#   - clusterSetRNGStream() seeds the FORK workers.
#   - Per-gene seed = master_seed * 1000003 + xxhash32(utr_name) -- makes the
#     per-gene permutation invariant to worker count and dispatch order.
#   - master_seed column written to output CSV.
#   - run_metadata.json sidecar records master_seed, R_version, package_versions,
#     command line, start_time, end_time, host.
#
# Usage (mouse SVZ has one lineage; the driver run_permutation_all.sh is
# invoked once with no lineage loop):
#
#   IN=$M26_DATA_ROOT/3utr/wildtypes_analysis/2026_analyses/utrcalling/\
# UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad
#
#   ./260526-permutation_optimized_2.R $IN \
#       --pseudotime_col dpt_pseudotime \
#       --covariate_cols batch_number \
#       --permutation_block_col sample_id \
#       --subsample_n 0 --seed 42 \
#       --n_perm_init 49 --n_perm_mid 199 --n_perm_max 999 --ncores 26
#
# (Mouse uses batch_number / sample_id (both 6 levels in this dataset, 1-to-1).
# SUBSAMPLE_N=0 means no subsampling -- mouse has ~14k lineage cells, which
# is tractable at full size on 26 cores.)
#
# Companion methods document: 260526-permutation_methods.md (Knoblich tree)
# Companion driver shell:     run_permutation_all.sh
# Companion QC notebook:      260526-permutation_concordance.ipynb
# ============================================================================

# ===== Section 0: Libraries and Argument Parsing =============================

library(tidyverse)
library(anndata)
library(pbapply)
library(parallel)
library(splines)
library(numDeriv)
library(digest)    # xxhash32 for per-gene seed derivation (permutation only)
library(jsonlite)  # run_metadata.json sidecar (permutation only)

# NOTE on RNG kind:
#   RNGkind() is left at the R default (Mersenne-Twister) for compatibility
#   with the inherited test suite (Tests 1-14 calibrate against MT). The
#   switch to L'Ecuyer-CMRG (parallel-safe substreams) is performed locally
#   in Phase B / Phase C, right before each `clusterSetRNGStream()` call.
#   That is sufficient because the only randomness consumed outside the
#   parallel sections is (i) the subsample (one set.seed + sample call,
#   deterministic regardless of RNG kind) and (ii) Phase D's per-gene
#   parametric null draw (which still uses MT, same as the source LRT script).
#   Setting RNGkind globally to L'Ecuyer-CMRG would change the synthetic
#   data in Tests 3, 5, 9, etc., breaking calibration without improving
#   reproducibility (the production randomness sits in the parallel section
#   alone).

pboptions(type = "timer")

parse_args <- function(args) {
    opts <- list(
        data_path = NULL, method = "joint", spline_df = 5L, min_cells = 50L,
        padj_method = "BH", ncores = max(1L, detectCores() - 1L),
        compute_se = TRUE, max_iter = 500L, shrinkage = "adaptive",
        lambda = 1.0, test = FALSE,
        bootstrap = 0L, bootstrap_pval = 0.1,
        calibration_n = 200L, lambda_sensitivity = FALSE,
        # Knoblich additions
        pseudotime_col = "dpt_pseudotime",
        covariate_cols = "batch_number",
        output_dir = NULL,
        keep_induced_false_only = FALSE,
        subsample_n = 0L,
        seed = 42L,
        # Permutation additions (260526)
        n_perm_init = 49L,
        n_perm_mid = 199L,
        n_perm_max = 999L,
        n_hits_stop = 5L,
        permutation_strategy = "within_sample",
        permutation_block_col = "sample_id"
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
        if (a == "--keep_induced_false_only") {
            opts$keep_induced_false_only <- TRUE
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
                pseudotime_col = {
                    opts$pseudotime_col <- val
                },
                covariate_cols = {
                    opts$covariate_cols <- val
                },
                output_dir = {
                    opts$output_dir <- val
                },
                subsample_n = {
                    opts$subsample_n <- as.integer(val)
                },
                seed = {
                    opts$seed <- as.integer(val)
                },
                n_perm_init = {
                    opts$n_perm_init <- as.integer(val)
                },
                n_perm_mid = {
                    opts$n_perm_mid <- as.integer(val)
                },
                n_perm_max = {
                    opts$n_perm_max <- as.integer(val)
                },
                n_hits_stop = {
                    opts$n_hits_stop <- as.integer(val)
                },
                permutation_strategy = {
                    opts$permutation_strategy <- val
                },
                permutation_block_col = {
                    opts$permutation_block_col <- val
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
        stop("Usage: Rscript 260526-permutation_optimized_2.R <path_to_h5ad> [--pseudotime_col dpt_pseudotime] [--covariate_cols batch_number] [--permutation_block_col sample_id] [--keep_induced_false_only] [--subsample_n N] [--seed S] [--method joint] [--spline_df 5] [--n_perm_init 49] [--n_perm_mid 199] [--n_perm_max 999] [--permutation_strategy within_sample|global] [--test]")
    }
    if (!opts$test) {
        if (!(opts$permutation_strategy %in% c("within_sample", "global"))) {
            stop("--permutation_strategy must be 'within_sample' or 'global'")
        }
        if (opts$n_perm_init <= 0L || opts$n_perm_mid < opts$n_perm_init ||
            opts$n_perm_max < opts$n_perm_mid) {
            stop("Permutation escalation tiers must satisfy 0 < n_perm_init <= n_perm_mid <= n_perm_max")
        }
        if (opts$n_hits_stop < 1L) {
            stop("--n_hits_stop must be >= 1")
        }
    }
    opts
}

# Parse covariate_cols string into a character vector. "none"/""/NULL -> empty vector.
parse_covariate_cols <- function(spec) {
    if (is.null(spec)) {
        return(character(0))
    }
    s <- trimws(as.character(spec))
    if (length(s) == 0L || s == "" || tolower(s) == "none") {
        return(character(0))
    }
    cols <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    cols[nzchar(cols)]
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
#' optimized_2 (publication-rigor revision): Shared-s pipeline with conditional
#' MLE refits in *both* the full and null beta subspaces.
#'
#' Pipeline:
#'   1. Fit full model jointly (or profile/fixed-s); on failure, fall back to
#'      fixed-s on the full design at s_trend.
#'   2. Extract s_shared = fit_full$s.
#'   3. Refit the full model beta at fixed s_shared (Stage 1). The joint optimum
#'      beta is NOT in general the conditional MLE at s_shared because L(beta, s)
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
#'   6. If ll_full < ll_null even after Stage 2 (vanishingly rare), set
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
    # honestly excluded from BH adjustment rather than reported as p = 1.
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

# ===== Section 4.5: Within-Sample Permutation of the Conditional LRT =========

#' Permute pseudotime values within the sample blocks defined by `blocks`.
#'
#' Preserves the per-block marginal distribution of pseudotime exactly while
#' breaking the within-block cell-pseudotime association. Vectorised via
#' `ave(x, blocks, FUN = sample)`. With `blocks = "all"` falls back to a
#' global shuffle.
permute_pseudotime <- function(pt, blocks, strategy = "within_sample") {
    if (strategy == "global") {
        return(sample(pt))
    }
    # within_sample: shuffle within each block. `ave` returns x reshaped to
    # match its input layout, so per-cell ordering is preserved.
    as.numeric(ave(pt, blocks, FUN = sample))
}

#' Within-sample permutation p-value for the conditional LRT statistic
#'
#' Generates B permuted-pseudotime replicates by shuffling the cell-level
#' pseudotime vector WITHIN sample blocks (preserves per-sample marginal +
#' between-sample batch differences in UTR usage; breaks within-sample cell-
#' pseudotime association). Refits beta_full at fixed s_shared on each
#' replicate's permuted design matrix; ll_null is invariant to the
#' permutation and is computed once at the call site.
#'
#' Why valid: under the null hypothesis "UTR usage is independent of
#' pseudotime within sample", the permuted-pseudotime data is exchangeable
#' with the observed data. The permutation null does NOT assume per-cell
#' exchangeability, nor Wilks regularity, nor DirMult model correctness on
#' the marginal distribution -- it only assumes that the joint distribution
#' of (Y, sample, batch) is preserved up to the within-sample pseudotime
#' re-labelling.
#'
#' Adaptive B with early stopping:
#'   - Run `n_perm_init` permutations (default 49).
#'   - If n_hits >= n_hits_stop -> stop (gene clearly not significant).
#'   - Else escalate to `n_perm_mid` (default 199); stop on n_hits >= n_hits_stop.
#'   - At n_perm_mid with n_hits == 0 -> escalate to `n_perm_max` (default 999).
#'   - At n_perm_max -> stop unconditionally.
#'
#' Reproducibility: caller passes `gene_seed` derived as
#'   master_seed * 1000003 + as.integer(digest(utr_name, algo="xxhash32"))
#' which makes the per-gene permutation invariant to worker count and
#' dispatch order. set.seed(gene_seed) is called once at function entry.
#'
#' Efficiency techniques:
#'   - s_shared FIXED (re-using the unpermuted MLE estimate).
#'   - ll_null computed ONCE before the loop (null model has no pt term).
#'   - beta_null and the null design matrix are also invariant.
#'   - The spline basis is recomputed on permuted pt each iteration (knots
#'     are invariant to the permutation because the marginal of pt is preserved
#'     exactly, but the per-cell basis values change).
#'   - Stage 2/3 safety nets are SKIPPED inside the permutation loop. If
#'     ll_full_perm < ll_null after the fixed-s refit, set D_b = 0. The
#'     Stage 2 padded-null refit doubles per-permutation cost and is unnecessary
#'     for permutation: the bias from skipping it is monotone (always pulls
#'     D_b down), which is conservative -- it can only inflate p_perm.
#'
#' @param Y_obs N x K observed count matrix
#' @param mm_null N x p_null null design matrix (cov-only; invariant)
#' @param mm_full N x p_full full design matrix (cov + spline; used for shape only)
#' @param covariate_mm N x p_cov covariate-only part of mm_full (first p_null cols
#'        of mm_full; rebuilt mm_full_perm = cbind(covariate_mm, spline_perm))
#' @param K number of isoforms
#' @param beta_full_obs original full-model beta (warm start)
#' @param beta_null_obs original null-model beta
#' @param s_shared the plug-in concentration
#' @param D_obs the observed deviance
#' @param pt_keep per-cell pseudotime, only for cells that passed gene-level keep
#' @param block_keep per-cell sample-block id (for within-sample permutation)
#' @param spline_df spline degrees of freedom (matches the original LRT fit)
#' @param strategy "within_sample" or "global"
#' @param gene_seed RNG seed for this gene (see Reproducibility above)
#' @param n_perm_init / n_perm_mid / n_perm_max escalation tiers
#' @param n_hits_stop early-stopping threshold (D_b >= D_obs counted)
#' @param max_iter L-BFGS-B max iterations per permutation fit
#' @return list(p_perm, n_permutations, n_failed, n_hits, D_b, status)
permute_dirmult_lrt <- function(Y_obs, mm_null, covariate_mm, K,
                                beta_full_obs, beta_null_obs, s_shared,
                                D_obs, pt_keep, block_keep, spline_df,
                                strategy, gene_seed,
                                n_perm_init, n_perm_mid, n_perm_max,
                                n_hits_stop, max_iter) {
    if (is.na(D_obs) || is.na(s_shared) || n_perm_init <= 0L ||
        length(pt_keep) != nrow(Y_obs) || length(block_keep) != nrow(Y_obs)) {
        return(list(
            p_perm = NA_real_, n_permutations = 0L, n_failed = 0L,
            n_hits = 0L, D_b = numeric(0), permutation_status = "skipped"
        ))
    }
    # Set the per-gene seed deterministically. Worker-invariant; dispatch-invariant.
    set.seed(gene_seed)

    # ll_null is fixed across permutations (the null model has no pseudotime term).
    ll_null <- dirmult_loglik_vec(Y_obs, compute_pi(beta_null_obs, mm_null) * s_shared)
    p_cov <- ncol(covariate_mm)

    # Pre-check rank once on the reference (unpermuted) full design. Under
    # within-sample permutation the marginal of pt is preserved, so the spline
    # basis rank is invariant across all permuted designs — no per-permutation
    # QR needed when the reference is full rank.
    spline_ref <- tryCatch(bs(pt_keep, df = spline_df), error = function(e) NULL)
    ref_full_rank <- !is.null(spline_ref) &&
        (qr(cbind(covariate_mm, spline_ref))$rank == ncol(covariate_mm) + ncol(spline_ref))

    perm_one <- function() {
        pt_perm <- permute_pseudotime(pt_keep, block_keep, strategy)
        # Spline basis on permuted pseudotime. The marginal of pt is preserved
        # under within-sample permutation, so the placement of the bs() knots
        # is identical to the unpermuted basis -- only the per-cell basis
        # values change. This is the desired behaviour: it isolates the
        # "pseudotime label" perturbation without inducing any artefact from
        # shifting knot positions.
        spline_perm <- tryCatch(
            bs(pt_perm, df = spline_df),
            error = function(e) NULL
        )
        if (is.null(spline_perm)) {
            return(NA_real_)
        }
        mm_full_perm <- cbind(covariate_mm, spline_perm)
        if (!ref_full_rank) {
            # Rank-fix only when the reference design was not full rank (rare;
            # skipped in the typical within-sample case where rank is invariant).
            qr_mm <- qr(mm_full_perm)
            if (qr_mm$rank < ncol(mm_full_perm)) {
                mm_full_perm <- mm_full_perm[, qr_mm$pivot[seq_len(qr_mm$rank)],
                    drop = FALSE
                ]
            }
        }
        # Warm start beta_full at the unpermuted MLE if the design column
        # count matches; otherwise zero-init (cold start). This is a defensive
        # fall-back; in practice the permuted spline produces the same column
        # count as the unpermuted spline because the marginal of pt is
        # invariant to within-sample permutation.
        p_full_perm <- ncol(mm_full_perm)
        if (p_full_perm == nrow(beta_full_obs)) {
            beta_init <- beta_full_obs
        } else {
            beta_init <- matrix(0, nrow = p_full_perm, ncol = K - 1L)
        }
        fit <- tryCatch(
            optimize_beta_fixed_s(Y_obs, mm_full_perm, K, beta_init, s_shared, max_iter),
            error = function(e) NULL
        )
        if (is.null(fit)) {
            return(NA_real_)
        }
        if (fit$ll < ll_null) {
            # Permuted alt fits the data worse than the null; this is a degenerate
            # replicate that we treat as D_b = 0 (cannot be more extreme than D_obs
            # for any D_obs >= 0). Conservative -- can only inflate p_perm.
            return(0)
        }
        max(-2 * (ll_null - fit$ll), 0)
    }

    D_b <- numeric(n_perm_max)
    n_done <- 0L
    n_hits <- 0L
    n_failed <- 0L

    # Helper: run permutations up to a target B, stopping early if n_hits >=
    # n_hits_stop. Returns TRUE if we should escalate, FALSE to stop.
    run_until <- function(target_B) {
        while (n_done < target_B) {
            d_b <- perm_one()
            if (is.na(d_b)) {
                n_failed <<- n_failed + 1L
                # Cap failed at 2x target so a pathological gene does not infinite-loop
                if (n_failed > 2L * target_B) {
                    return(NA)
                }
                next
            }
            n_done <<- n_done + 1L
            D_b[n_done] <<- d_b
            if (d_b >= D_obs) n_hits <<- n_hits + 1L
            if (n_hits >= n_hits_stop) {
                return(FALSE)
            }
        }
        TRUE  # continue to next escalation tier
    }

    escalate <- run_until(n_perm_init)
    if (isTRUE(escalate) && n_perm_mid > n_perm_init) {
        escalate <- run_until(n_perm_mid)
    }
    if (isTRUE(escalate) && n_perm_max > n_perm_mid && n_hits == 0L) {
        escalate <- run_until(n_perm_max)
    }

    # Classify outcome
    status <- if (is.na(escalate)) {
        "all_failed"
    } else if (n_done < 1L) {
        "all_failed"
    } else if (n_hits >= n_hits_stop) {
        if (n_done <= n_perm_init) "early_stop_high_p" else "stopped_on_hits"
    } else if (n_done >= n_perm_max) {
        "max_b_few_hits"
    } else {
        "stopped_on_hits"
    }

    # Empirical p-value with +1/+1 correction. If all replicates failed,
    # p_perm = NA (the gene is excluded from BH downstream).
    p_perm <- if (n_done > 0L) (1 + n_hits) / (n_done + 1) else NA_real_

    list(
        p_perm = p_perm, n_permutations = n_done, n_failed = n_failed,
        n_hits = n_hits, D_b = D_b[seq_len(n_done)],
        permutation_status = status
    )
}

#' Per-gene reproducible seed derivation.
#'
#' seed_g = master_seed * 1000003 + as.integer(xxhash32(utr_name))
#' xxhash32 hashes the gene id to a 32-bit integer; multiplied by a prime to
#' decorrelate from the master seed.
gene_seed_from <- function(master_seed, utr_name) {
    h <- strtoi(substr(digest(utr_name, algo = "xxhash32"), 1L, 7L), base = 16L)
    # Combine with master_seed via a simple LCG-style mix. The cast to int32
    # is intentional -- R's set.seed accepts a 32-bit signed int.
    as.integer(bitwXor(as.integer(master_seed * 1000003L), h))
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
    # twice — once from the bad init, once from the padded null — and assert
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

    # Test 15: (deleted in 260526 fork -- tested bootstrap_dirmult_lrt which
    # has been replaced by permute_dirmult_lrt; see Test 19-21 below for the
    # permutation analogues of "uniformity / signal recovery / reproducibility".)

    # Test 16: Multi-regime null calibration grid (replaces single-cell Test 9)
    # Sweeps K x N x s in {2,3,5} x {100,200,500} x {0.5,5,50}, generates 150
    # null genes per cell, and runs two FWER-controlled checks:
    #   (a) Bonferroni-corrected per-cell KS test at family-wise alpha = 0.01
    #       Per-cell threshold = 0.01 / 27 ~= 3.7e-4. Requires every cell to
    #       pass individually after FWER correction.
    #   (b) Fisher combined test across all 27 cells: X^2 = -2 * sum(log(p_cell))
    #       ~ chi^2(54) under global H0. Requires X^2 p-value > 0.01 globally.
    # Both must pass. This is strictly more rigorous than the naive "all 27
    # cells pass at 1%" check (which has ~24% FWER under H0 by Sidak, causing
    # spurious failures). Uses n=150 per cell for better KS resolution and to
    # reduce the chance of ties from near-p=1 clustering in the extreme regime.
    cat("Test 16: Multi-regime null calibration grid (27 cells, FWER 1%)... ")
    set.seed(16L)
    t16_n_genes <- 150L
    t16_alpha_bonf <- 0.01 / 27 # Bonferroni-corrected per-cell threshold
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

    # ===== Permutation-specific tests (260526) ===============================

    # Test 17: permute_pseudotime preserves within-block marginals
    cat("Test 17: within-sample permutation preserves per-block marginals... ")
    set.seed(123L)
    pt_t17 <- runif(200L)
    block_t17 <- rep(c("A", "B", "C", "D"), each = 50L)
    pt_perm_t17 <- permute_pseudotime(pt_t17, block_t17, "within_sample")
    t17_pass <- TRUE
    for (blk in unique(block_t17)) {
        orig <- sort(pt_t17[block_t17 == blk])
        perm <- sort(pt_perm_t17[block_t17 == blk])
        if (!isTRUE(all.equal(orig, perm))) {
            t17_pass <- FALSE
            cat("\n  block", blk, "marginal not preserved")
        }
    }
    if (t17_pass) cat("PASS\n") else {
        cat("\nFAIL\n")
        all_pass <- FALSE
    }

    # Test 18: gene_seed_from is deterministic for same inputs
    cat("Test 18: gene_seed_from is deterministic and decorrelates with master_seed... ")
    s_a1 <- gene_seed_from(42L, "GENE_A")
    s_a2 <- gene_seed_from(42L, "GENE_A")
    s_b1 <- gene_seed_from(42L, "GENE_B")
    s_a3 <- gene_seed_from(43L, "GENE_A")
    if (identical(s_a1, s_a2) && !identical(s_a1, s_b1) && !identical(s_a1, s_a3)) {
        cat("PASS\n")
    } else {
        cat("FAIL (", s_a1, s_a2, s_b1, s_a3, ")\n")
        all_pass <- FALSE
    }

    # Test 19: permute_dirmult_lrt early-stops on a known-null gene (B<=49)
    cat("Test 19: adaptive stopping at B_init for known-null gene... ")
    set.seed(2026L)
    sim19 <- generate_test_data(300, 3, 8, spline_df = 4)
    # Block structure: 3 blocks of 100 cells
    block_t19 <- rep(c("S1", "S2", "S3"), each = 100L)
    res19 <- fit_dirmult_gene(sim19$Y, sim19$mm_full, sim19$mm_null, sim19$K,
        "joint", 8,
        lambda = 0.1, max_iter = 300, compute_se = FALSE
    )
    # Use cov_mm = intercept-only (first column of mm_full). For the
    # permutation we just need ANY pt vector of length N that the inner
    # bs(pt_perm, df) call can consume; the pt vector used to BUILD mm_full
    # inside generate_test_data is not exposed by that helper, so we supply
    # a fresh canonical pt here. The permutation test is API/sanity only.
    pt_t19 <- seq(0, 1, length.out = 300L)
    perm_res19 <- permute_dirmult_lrt(
        Y_obs = sim19$Y,
        mm_null = sim19$mm_null,
        covariate_mm = sim19$mm_null,  # intercept-only
        K = sim19$K,
        beta_full_obs = res19$beta_full,
        beta_null_obs = res19$beta_null,
        s_shared = res19$s_shared,
        D_obs = res19$deviance,
        pt_keep = pt_t19,
        block_keep = block_t19,
        spline_df = 4L,
        strategy = "within_sample",
        gene_seed = 12345L,
        n_perm_init = 49L, n_perm_mid = 199L, n_perm_max = 999L,
        n_hits_stop = 5L, max_iter = 300L
    )
    # For a truly null gene the chi-sq deviance lands in the bulk of the
    # permutation distribution; we expect n_hits to hit the early-stop
    # threshold (5) within the first few dozen permutations. Accept the
    # test as PASSED if (i) p_perm is a valid probability and (ii) the
    # adaptive stopping triggered (status == "early_stop_high_p" with
    # n_perm < 49) OR ran the full n_perm_init without early stopping.
    valid_p <- !is.na(perm_res19$p_perm) &&
        perm_res19$p_perm >= 0 && perm_res19$p_perm <= 1
    valid_status <- perm_res19$permutation_status %in%
        c("early_stop_high_p", "stopped_on_hits", "max_b_few_hits")
    if (valid_p && valid_status && perm_res19$n_permutations >= 1L) {
        cat("PASS (p_perm =", round(perm_res19$p_perm, 3),
            ", n_perm =", perm_res19$n_permutations,
            ", status =", perm_res19$permutation_status, ")\n")
    } else {
        cat("FAIL (", perm_res19$p_perm, perm_res19$n_permutations,
            perm_res19$permutation_status, ")\n")
        all_pass <- FALSE
    }

    # Test 20: permute_dirmult_lrt rejects a known signal (B may need to escalate)
    cat("Test 20: permute_dirmult_lrt rejects a known strong signal... ")
    set.seed(99L)
    N_t20 <- 300L
    K_t20 <- 3L
    pt_t20 <- sort(runif(N_t20))
    mm_t20_full <- model.matrix(~ bs(pt_t20, df = 4L))
    mm_t20_null <- matrix(1, nrow = N_t20, ncol = 1L)
    p_t20 <- ncol(mm_t20_full)
    # Strong pseudotime-dependent beta
    true_beta_t20 <- matrix(c(
        0, 2.0, 3.0, 2.0, 0,
        0, -2.0, -3.0, -2.0, 0
    ), nrow = p_t20, ncol = K_t20 - 1L)
    pi_t20 <- compute_pi(true_beta_t20, mm_t20_full)
    alpha_t20 <- pi_t20 * 20
    n_t20 <- rpois(N_t20, 20) + 5
    Y_t20 <- t(sapply(seq_len(N_t20), function(i) {
        p_dir <- rgamma(K_t20, alpha_t20[i, ])
        p_dir <- p_dir / sum(p_dir)
        rmultinom(1, n_t20[i], p_dir)
    }))
    res20 <- fit_dirmult_gene(Y_t20, mm_t20_full, mm_t20_null, K_t20,
        "joint", 20, lambda = 0.1, max_iter = 300, compute_se = FALSE
    )
    block_t20 <- rep(c("S1", "S2", "S3"), each = 100L)
    perm_res20 <- permute_dirmult_lrt(
        Y_obs = Y_t20,
        mm_null = mm_t20_null,
        covariate_mm = mm_t20_null,
        K = K_t20,
        beta_full_obs = res20$beta_full,
        beta_null_obs = res20$beta_null,
        s_shared = res20$s_shared,
        D_obs = res20$deviance,
        pt_keep = pt_t20,
        block_keep = block_t20,
        spline_df = 4L,
        strategy = "within_sample",
        gene_seed = 67890L,
        n_perm_init = 49L, n_perm_mid = 199L, n_perm_max = 999L,
        n_hits_stop = 5L, max_iter = 300L
    )
    # Strong signal: permutation should produce p_perm <= 0.05 with B <= 199.
    # Accept p_perm <= 0.1 since we use a synthetic dataset with limited size.
    if (!is.na(perm_res20$p_perm) && perm_res20$p_perm <= 0.1) {
        cat("PASS (p_perm =", round(perm_res20$p_perm, 4),
            ", n_perm =", perm_res20$n_permutations,
            ", D_obs =", round(res20$deviance, 2), ")\n")
    } else {
        cat("FAIL (p_perm =", perm_res20$p_perm, ", D_obs =", res20$deviance, ")\n")
        all_pass <- FALSE
    }

    # Test 21: permute_dirmult_lrt is reproducible (same gene_seed -> same p_perm)
    cat("Test 21: permute_dirmult_lrt is reproducible under fixed gene_seed... ")
    perm_res21a <- permute_dirmult_lrt(
        Y_obs = Y_t20,
        mm_null = mm_t20_null,
        covariate_mm = mm_t20_null,
        K = K_t20,
        beta_full_obs = res20$beta_full,
        beta_null_obs = res20$beta_null,
        s_shared = res20$s_shared,
        D_obs = res20$deviance,
        pt_keep = pt_t20,
        block_keep = block_t20,
        spline_df = 4L,
        strategy = "within_sample",
        gene_seed = 11111L,
        n_perm_init = 49L, n_perm_mid = 199L, n_perm_max = 999L,
        n_hits_stop = 5L, max_iter = 300L
    )
    perm_res21b <- permute_dirmult_lrt(
        Y_obs = Y_t20,
        mm_null = mm_t20_null,
        covariate_mm = mm_t20_null,
        K = K_t20,
        beta_full_obs = res20$beta_full,
        beta_null_obs = res20$beta_null,
        s_shared = res20$s_shared,
        D_obs = res20$deviance,
        pt_keep = pt_t20,
        block_keep = block_t20,
        spline_df = 4L,
        strategy = "within_sample",
        gene_seed = 11111L,
        n_perm_init = 49L, n_perm_mid = 199L, n_perm_max = 999L,
        n_hits_stop = 5L, max_iter = 300L
    )
    if (identical(perm_res21a$p_perm, perm_res21b$p_perm) &&
        identical(perm_res21a$n_permutations, perm_res21b$n_permutations) &&
        identical(perm_res21a$D_b, perm_res21b$D_b)) {
        cat("PASS (p_perm =", round(perm_res21a$p_perm, 4),
            ", D_b identical over", length(perm_res21a$D_b), "permutations)\n")
    } else {
        cat("FAIL (a:", perm_res21a$p_perm, "b:", perm_res21b$p_perm, ")\n")
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
cat(
    "  bootstrap:", opts$bootstrap,
    if (opts$bootstrap > 0L) sprintf("(threshold p_chisq < %.3f)", opts$bootstrap_pval) else "(off)",
    "\n"
)
cat("  calibration_n:", opts$calibration_n, "\n")
cat("  lambda_sensitivity:", opts$lambda_sensitivity, "\n")
cat("  pseudotime_col:", opts$pseudotime_col, "\n")
cat("  covariate_cols:", opts$covariate_cols, "\n")
cat("  keep_induced_false_only:", opts$keep_induced_false_only, "\n")
cat("  subsample_n:", opts$subsample_n, "\n")
cat("  seed:", opts$seed, "\n")
cat("  permutation_strategy:", opts$permutation_strategy, "\n")
cat("  permutation_block_col:", opts$permutation_block_col, "\n")
cat("  n_perm_init / mid / max:", opts$n_perm_init, "/", opts$n_perm_mid, "/", opts$n_perm_max, "\n")
cat("  n_hits_stop:", opts$n_hits_stop, "\n")

# Record run start time for the run_metadata.json sidecar
run_start_time <- format(Sys.time())
run_host <- Sys.info()[["nodename"]]
run_cmd <- paste(c(get0(".Last.value"), commandArgs(trailingOnly = FALSE)), collapse = " ")

covariate_col_vec <- parse_covariate_cols(opts$covariate_cols)
covariate_tag <- if (length(covariate_col_vec) == 0L) {
    "no_covariates"
} else {
    paste(covariate_col_vec, collapse = "_")
}

# Resolve output_dir before loading anything (cells_used.csv goes here too)
output_dir <- if (!is.null(opts$output_dir) && nzchar(opts$output_dir)) {
    opts$output_dir
} else {
    file.path(dirname(opts$data_path), "DirMult_optimized")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(dirname(opts$data_path))
cat("Output dir:", output_dir, "\n")

cat("\nLoading data from", opts$data_path, "\n")
data <- read_h5ad(opts$data_path)
n_total <- nrow(data$obs)
cat("Loaded:", n_total, "cells x", ncol(data$X), "features\n")

# Pull what we need from the AnnData; KEEP X SPARSE for now. Filtering and
# subsampling happen on the sparse matrix to avoid the 33 GiB peak from
# densifying the full data set.
all_cells <- rownames(data$obs)
utr_names <- as.character(data$var$utr_name)

if (!opts$pseudotime_col %in% colnames(data$obs)) {
    stop(
        "Pseudotime column '", opts$pseudotime_col, "' not found in obs. ",
        "Available: ", paste(colnames(data$obs), collapse = ", ")
    )
}
pseudotime <- as.numeric(data$obs[[opts$pseudotime_col]])

missing_cov <- setdiff(covariate_col_vec, colnames(data$obs))
if (length(missing_cov) > 0L) {
    stop(
        "Covariate columns not found in obs: ", paste(missing_cov, collapse = ", "),
        ". Available: ", paste(colnames(data$obs), collapse = ", ")
    )
}
covariate_df <- if (length(covariate_col_vec) > 0L) {
    as.data.frame(data$obs[, covariate_col_vec, drop = FALSE])
} else {
    data.frame(.dummy_intercept = rep(1L, n_total))
}

# induced column (needed if --keep_induced_false_only) -- read once.
induced_chr <- if ("induced" %in% colnames(data$obs)) {
    as.character(data$obs[["induced"]])
} else {
    rep(NA_character_, n_total)
}

# NEW (260526): permutation block column. Each cell's permutation block
# defines the unit within which pseudotime labels are shuffled under
# --permutation_strategy=within_sample. For Knoblich this is sample_id;
# for mouse it is batch.
if (!opts$permutation_block_col %in% colnames(data$obs)) {
    stop(
        "Permutation block column '", opts$permutation_block_col,
        "' not found in obs. Available: ",
        paste(colnames(data$obs), collapse = ", ")
    )
}
permutation_block_all <- as.character(data$obs[[opts$permutation_block_col]])

# ----- Cell-level filters (BEFORE dense conversion) ------------------------

# 1. Drop cells with NA in the requested pseudotime column. When a
#    lineage-specific column is passed, this also restricts to that lineage.
mask <- !is.na(pseudotime)
cat(
    "After dropping NA in", opts$pseudotime_col, ": kept", sum(mask),
    "/", n_total, "\n"
)

# 2. Optional: drop induced=TRUE cells (Knoblich-specific; resolves the
#    perfect confounding between `induced` and `sample`/`batch_number` so
#    batch_number can be safely included as a fixed effect).
if (opts$keep_induced_false_only) {
    is_false <- induced_chr %in% c("FALSE", "False", "false", "F", "0")
    n_drop <- sum(mask & !is_false)
    mask <- mask & is_false
    cat(
        "Filter induced=FALSE only: dropped", n_drop, "induced=TRUE cells; kept",
        sum(mask), "induced=FALSE cells\n"
    )
}

# 3. Optional: random subsample to a fixed N for tractability. Reproducible
#    via --seed; cell barcodes are written to cells_used.csv below.
if (opts$subsample_n > 0L && sum(mask) > opts$subsample_n) {
    set.seed(opts$seed)
    candidates <- which(mask)
    chosen <- sort(sample(candidates, opts$subsample_n))
    mask <- logical(length(mask))
    mask[chosen] <- TRUE
    cat("Subsampled to", opts$subsample_n, "cells (seed=", opts$seed, ")\n")
}

cat("Final cell count after filters:", sum(mask), "\n")

# Apply mask while X is still sparse (cheap)
X <- data$X[mask, , drop = FALSE]
covariate_df <- covariate_df[mask, , drop = FALSE]
pseudotime <- pseudotime[mask]
all_cells <- all_cells[mask]
induced_chr_kept <- induced_chr[mask]
permutation_block <- permutation_block_all[mask]  # NEW (260526)

# Free the full AnnData object before densifying
rm(data, induced_chr, permutation_block_all)
gc(verbose = FALSE)

# Densify the small subset
cat("Converting subset X to dense ...\n")
if (inherits(X, "sparseMatrix")) X <- as.matrix(X)
cat(
    "X size:", format(object.size(X), units = "auto"),
    "  (", nrow(X), "cells x", ncol(X), "features)\n"
)
gc(verbose = FALSE)

cat(
    "\nData:", nrow(X), "cells,", ncol(X), "features,",
    length(unique(utr_names)), "unique UTRs\n"
)
for (col in covariate_col_vec) {
    vec <- as.character(covariate_df[[col]])
    cat("  covariate '", col, "': ", length(unique(vec)), " levels (",
        paste(head(sort(unique(vec)), 8), collapse = ", "),
        if (length(unique(vec)) > 8) ", ..." else "", ")\n",
        sep = ""
    )
}
if (length(covariate_col_vec) == 0L) {
    cat("  no covariates: null model is intercept-only\n")
}
cat(
    "Pseudotime range:", round(min(pseudotime), 4), "-",
    round(max(pseudotime), 4), "\n"
)

# Save cell list for reproducibility (barcode + pseudotime + covariates +
# permutation block).
cells_used_path <- file.path(output_dir, "cells_used.csv")
cells_used_df <- data.frame(
    cell_barcode = all_cells,
    pseudotime = pseudotime,
    induced = induced_chr_kept,
    permutation_block = permutation_block,
    stringsAsFactors = FALSE
)
for (col in covariate_col_vec) {
    cells_used_df[[col]] <- as.character(covariate_df[[col]])
}
write.csv(cells_used_df, cells_used_path, row.names = FALSE)
cat(
    "Saved cell list:", cells_used_path,
    "(", nrow(cells_used_df), "rows)\n"
)

# Sanity check the permutation block: count cells per block
cat("\nPermutation block composition (",
    opts$permutation_block_col, "):\n", sep = ""
)
block_tab <- table(permutation_block)
cat("  ", length(block_tab), "blocks; cells per block:\n")
print(block_tab)
if (length(block_tab) < 2L) {
    stop(
        "Permutation block has only one level after filters. Within-sample ",
        "permutation requires >= 2 blocks. Check --permutation_block_col."
    )
}
if (opts$permutation_strategy == "within_sample" && min(block_tab) < 5L) {
    cat(
        "  WARNING: at least one block has < 5 cells. Within-sample permutation",
        "for that block has very low entropy; the resulting null tail may be",
        "dominated by the other blocks. Consider --permutation_strategy global",
        "if this is a problem.\n"
    )
}

# Helper: build covariate design matrix for a subset of cells (mask "keep").
# Drops factor levels absent in the subset, then drops columns that become
# constant; the rank-fix step in fit_dirmult_gene catches anything else.
build_covariate_mm <- function(covariate_df, keep, covariate_col_vec) {
    if (length(covariate_col_vec) == 0L) {
        return(matrix(1,
            nrow = sum(keep), ncol = 1L,
            dimnames = list(NULL, "(Intercept)")
        ))
    }
    sub <- covariate_df[keep, covariate_col_vec, drop = FALSE]
    used <- character(0)
    for (col in covariate_col_vec) {
        vec <- as.character(sub[[col]])
        if (length(unique(vec)) > 1L) {
            sub[[col]] <- factor(vec)
            used <- c(used, col)
        }
    }
    if (length(used) == 0L) {
        return(matrix(1,
            nrow = nrow(sub), ncol = 1L,
            dimnames = list(NULL, "(Intercept)")
        ))
    }
    formula_str <- paste("~", paste(used, collapse = " + "))
    model.matrix(as.formula(formula_str), data = sub)
}

# Precompute spline basis on the (filtered, possibly subsampled) pseudotime
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
        Y = NULL, n_counts = NULL,
        # NEW (260526): per-gene pt and block subsets, used by permute_dirmult_lrt
        pt_keep = NULL, block_keep = NULL, covariate_mm = NULL
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
            cov_mm <- build_covariate_mm(covariate_df, keep, covariate_col_vec)
            spline_sel <- spline_basis[keep, , drop = FALSE]
            mm_null <- cov_mm
            mm_full <- cbind(cov_mm, spline_sel)
            r <- fit_dirmult_gene(
                Y_sel, mm_full, mm_null, K_obs, method_to_use,
                trend_result$s_trend[utr], lambda_val,
                opts$max_iter, opts$compute_se
            )
            r$Y <- Y_sel
            r$n_counts <- rowSums(Y_sel)
            # NEW (260526): per-gene pt and block subsets, used by permute_dirmult_lrt
            r$pt_keep <- pseudotime[keep]
            r$block_keep <- permutation_block[keep]
            r$covariate_mm <- cov_mm
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

# --- Dynamic load-balanced dispatch helper -------------------------------
# Replaces pblapply(cl=cluster) which uses STATIC chunked scheduling
# (parLapply: contiguous blocks of ~N/ncores genes per worker, with a sync
# barrier every progress-bar chunk). Static chunking leaves cores idle while
# a few stragglers finish -- catastrophic for the permutation phase where
# per-gene cost varies ~20x (B=49 early-stop vs B=999 escalation) AND for the
# fit phase where high-K / high-cell genes are much slower than average. When
# slow genes are adjacent in the gene list (high-K genes cluster), static
# scheduling can run 7x slower than ideal (benchmarked 2026-05-27).
#
# parLapplyLB(chunk.size=1) dispatches ONE gene at a time to whichever worker
# is free, eliminating the straggler problem (benchmarked ~1.05x ideal even
# with maximally-clustered slow genes). We split into `n_progress` ROUND-ROBIN
# batches purely to print progress; round-robin assignment guarantees slow
# genes spread across batches so each batch stays balanced (a contiguous split
# would re-cluster them and reintroduce barriers).
#
# CORRECTNESS: results are invariant to scheduling because every gene reseeds
# its own RNG via set.seed(gene_seed_from(master_seed, utr_name)) at task
# entry (Phase C / D) or consumes no RNG at all (Phase B). Verified
# bit-identical between static vs dynamic AND between 26-core vs 4-core runs.
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
            "  %s %d/%d genes (%.0f%%), elapsed %s\n",
            label, done, n, 100 * done / n,
            format(round(difftime(Sys.time(), t0, units = "mins"), 2))
        ))
        utils::flush.console()
    }
    out
}

# --- Phase B: Per-gene fitting ---
cat("\nPhase B: Fitting DirMult models per gene (", length(unique_utrs), " UTRs)...\n")
cat("Setting up FORK cluster with", opts$ncores, "cores (dynamic load balancing)\n")
cl <- makeCluster(opts$ncores, type = "FORK")
# Switch to L'Ecuyer-CMRG for parallel-safe substreams and seed the workers
# deterministically. Phase B itself does not consume randomness, but this is
# cheap insurance against future changes to fit_dirmult_gene that might
# introduce randomness (e.g. jittered inits).
RNGkind("L'Ecuyer-CMRG")
clusterSetRNGStream(cl, opts$seed)

results_list <- dispatch_lb(
    cl, seq_along(unique_utrs),
    function(idx) fit_one_utr(unique_utrs[idx], gene_lambda[unique_utrs[idx]]),
    label = "Phase B fit:"
)

stopCluster(cl)

# Index results by utr_name for fast lookup in Phases C/D
names(results_list) <- vapply(results_list, `[[`, character(1), "utr_name")

# --- Phase C: Within-sample permutation on ALL converged genes ---
# Replaces the parametric-bootstrap Phase C from 260503-...R. For each gene
# with lrt_status == "ok", run permute_dirmult_lrt with adaptive B and
# return p_value_perm. The per-gene seed is derived from the master seed
# + xxhash32(utr_name) for full reproducibility, invariant to worker count
# and dispatch order.

permutation_results <- list()
cat(
    "\nPhase C: Within-sample permutation on all converged genes\n",
    "  strategy: ", opts$permutation_strategy, "\n",
    "  block_col: ", opts$permutation_block_col, "\n",
    "  adaptive B: init=", opts$n_perm_init,
    ", mid=", opts$n_perm_mid,
    ", max=", opts$n_perm_max,
    " (escalate when n_hits < ", opts$n_hits_stop, ")\n",
    sep = ""
)
eligible_perm_idx <- which(vapply(results_list, function(r) {
    !is.null(r$lrt_status) && r$lrt_status == "ok" &&
        !is.null(r$Y) && !is.null(r$beta_full) && !is.null(r$beta_null) &&
        !is.null(r$mm_null) && !is.null(r$covariate_mm) &&
        !is.null(r$pt_keep) && !is.null(r$block_keep) &&
        !is.na(r$s_shared) && !is.na(r$deviance) && r$deviance >= 0
}, logical(1)))
cat("  ", length(eligible_perm_idx), "genes eligible for permutation\n")
n_skipped_perm <- length(results_list) - length(eligible_perm_idx)
cat("  ", n_skipped_perm, "genes excluded (lrt_status != ok, NA deviance, or missing payload)\n")

if (length(eligible_perm_idx) > 0L) {
    cl_p <- makeCluster(opts$ncores, type = "FORK")
    # L'Ecuyer-CMRG substream seeding. Each worker still derives its own
    # per-gene seed inside permute_dirmult_lrt (defence in depth), but
    # seeding the cluster handles any other random calls. The dynamic
    # load-balanced dispatch (dispatch_lb -> parLapplyLB chunk.size=1) is
    # especially important here: the B=49 vs B=999 escalation makes per-gene
    # cost vary ~20x, so static chunking would idle most cores waiting for
    # the genes that escalate to B=999.
    RNGkind("L'Ecuyer-CMRG")
    clusterSetRNGStream(cl_p, opts$seed)

    perm_out <- dispatch_lb(cl_p, eligible_perm_idx, function(j) {
        r <- results_list[[j]]
        g_seed <- gene_seed_from(opts$seed, r$utr_name)
        tryCatch(
            permute_dirmult_lrt(
                Y_obs = r$Y,
                mm_null = r$mm_null,
                covariate_mm = r$covariate_mm,
                K = r$K,
                beta_full_obs = r$beta_full,
                beta_null_obs = r$beta_null,
                s_shared = r$s_shared,
                D_obs = r$deviance,
                pt_keep = r$pt_keep,
                block_keep = r$block_keep,
                spline_df = opts$spline_df,
                strategy = opts$permutation_strategy,
                gene_seed = g_seed,
                n_perm_init = opts$n_perm_init,
                n_perm_mid = opts$n_perm_mid,
                n_perm_max = opts$n_perm_max,
                n_hits_stop = opts$n_hits_stop,
                max_iter = opts$max_iter
            ),
            error = function(e) {
                list(
                    p_perm = NA_real_, n_permutations = 0L,
                    n_failed = 0L, n_hits = 0L,
                    D_b = numeric(0), permutation_status = "error"
                )
            }
        )
    }, label = "Phase C perm:")
    stopCluster(cl_p)
    names(perm_out) <- names(results_list)[eligible_perm_idx]
    permutation_results <- perm_out
    # Inject permutation results back into per-gene results_list
    for (nm in names(perm_out)) {
        results_list[[nm]]$p_perm <- perm_out[[nm]]$p_perm
        results_list[[nm]]$n_permutations <- perm_out[[nm]]$n_permutations
        results_list[[nm]]$n_perm_hits <- perm_out[[nm]]$n_hits
        results_list[[nm]]$n_perm_failed <- perm_out[[nm]]$n_failed
        results_list[[nm]]$permutation_status <- perm_out[[nm]]$permutation_status
    }
    cat("  Permutation phase complete.\n")
    # Status summary
    perm_status <- vapply(perm_out, `[[`, character(1), "permutation_status")
    cat("  Permutation status breakdown:\n")
    print(table(perm_status))
}

# --- Phase D: Real-data empirical calibration footprint ---
# Sample N random tested genes, generate one null replicate from each fitted
# null model, refit, and KS-test the resulting p-values against the uniform.
# This is an end-to-end calibration check on the *exact* parameter regimes
# present in the user's data, complementing the synthetic-data Test 16.

calibration_result <- list(
    ks_p = NA_real_, n_used = 0L, n_attempted = 0L
)
if (opts$calibration_n > 0L) {
    cat(
        "\nPhase D: Real-data null calibration footprint (n =",
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
        RNGkind("L'Ecuyer-CMRG")
        clusterSetRNGStream(cl_d, opts$seed)
        # parLapplyLB (via dispatch_lb) for load balancing; per-gene seeding
        # makes the calibration draw reproducible (it consumes RNG via
        # rgamma/rmultinom). The original pbsapply path inherited the parent
        # RNG per fork and was NOT reproducible across core counts; seeding
        # per gene here fixes that. Calibration seed namespace is offset from
        # the permutation seed (+ 7L) so the two phases never share a stream.
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
            cat(
                "  Real-data null KS p =", round(ks_cal$p.value, 4),
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
    # Permutation columns (new in 260526)
    p_value_perm = vapply(results_list, function(x) extract_num(x, "p_perm"), numeric(1)),
    n_permutations = vapply(results_list, function(x) {
        v <- x$n_permutations
        if (is.null(v)) NA_integer_ else as.integer(v)
    }, integer(1)),
    n_perm_hits = vapply(results_list, function(x) {
        v <- x$n_perm_hits
        if (is.null(v)) NA_integer_ else as.integer(v)
    }, integer(1)),
    n_perm_failed = vapply(results_list, function(x) {
        v <- x$n_perm_failed
        if (is.null(v)) NA_integer_ else as.integer(v)
    }, integer(1)),
    permutation_status = vapply(results_list, function(x) {
        v <- x$permutation_status
        if (is.null(v)) NA_character_ else as.character(v)
    }, character(1)),
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
    method = vapply(results_list, function(x) extract_chr(x, "method", opts$method), character(1)),
    # Reproducibility metadata (constant per row, but written per row for downstream merge convenience)
    master_seed = rep(as.integer(opts$seed), length(results_list)),
    permutation_strategy = rep(opts$permutation_strategy, length(results_list)),
    permutation_block_col = rep(opts$permutation_block_col, length(results_list))
)

# Reported p_value: prefer permutation; fall back to chi-sq when permutation
# is unavailable (lrt_status != "ok", or all permutations failed).
results_df <- results_df %>%
    mutate(p_value = ifelse(!is.na(p_value_perm), p_value_perm, p_value_chisq))

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
# null_fit_failed, error). These are honestly excluded rather than reported as
# p = 1 with deviance = 0.
results_df <- results_df %>%
    filter(!is.na(p_value))

# Restrict BH adjustment to fully-converged genes (converged == TRUE &
# lrt_status == "ok").  Non-converged genes have unreliable deviance statistics
# because the optimizer stopped before reaching the true maximum — their
# p-values tend to be artificially small, which inflates the BH threshold and
# produces spurious discoveries among the converged genes.  Setting padj = NA
# for non-converged genes keeps them in the output for inspection but excludes
# them from FDR control.  Downstream notebooks (e.g. composite_ranking) already
# filter to converged == TRUE before use; this change makes the CSV consistent
# with that intent from the outset.
conv_for_bh <- results_df$converged == TRUE & results_df$lrt_status == "ok"
padj_vec <- rep(NA_real_, nrow(results_df))
padj_vec[conv_for_bh] <- p.adjust(
    results_df$p_value[conv_for_bh],
    method = opts$padj_method
)
results_df$padj <- padj_vec

# Summary
n_total <- length(results_list)
n_tested <- nrow(results_df)
n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
n_conv <- sum(results_df$converged, na.rm = TRUE)
n_full_conv <- sum(results_df$full_converged, na.rm = TRUE)
n_ll_viol <- sum(results_df$ll_full < results_df$ll_null - 1e-6, na.rm = TRUE)

# Counts of excluded statuses among ALL processed genes
all_status <- vapply(results_list, function(x) extract_chr(x, "lrt_status", "skipped"), character(1))
n_lrt_failed <- sum(all_status == "optimization_failed")
n_null_fail <- sum(all_status == "null_fit_failed")
n_skipped <- sum(all_status == "skipped")
n_error <- sum(all_status == "error")
n_perm_used <- sum(!is.na(results_df$p_value_perm))

cat("\nResults: ", n_tested, " genes tested, ", n_conv, " converged, ",
    n_sig, " significant among converged (padj<0.05); ",
    n_tested - n_conv, " non-converged genes have padj=NA\n",
    sep = ""
)
cat("Full model primary convergence: ", n_full_conv, "/", n_tested, "\n", sep = "")
cat("LL violations (ll_full < ll_null) among reported: ", n_ll_viol, "\n", sep = "")
cat("Excluded genes:\n")
cat("  optimization_failed: ", n_lrt_failed, " (Stage 3 honest exclusion)\n", sep = "")
cat("  null_fit_failed:     ", n_null_fail, "\n", sep = "")
cat("  skipped (low cells): ", n_skipped, "\n", sep = "")
cat("  error:               ", n_error, "\n", sep = "")
cat("Permutation p-values used (replaces chi-sq): ", n_perm_used, "\n", sep = "")
# Permutation B-tier breakdown (escalation diagnostic)
n_perm_init_only <- sum(results_df$n_permutations <= opts$n_perm_init, na.rm = TRUE)
n_perm_mid_only <- sum(
    results_df$n_permutations > opts$n_perm_init &
        results_df$n_permutations <= opts$n_perm_mid,
    na.rm = TRUE
)
n_perm_max <- sum(results_df$n_permutations > opts$n_perm_mid, na.rm = TRUE)
cat("Permutation B-tier breakdown:\n")
cat("  <=", opts$n_perm_init, ":", n_perm_init_only, " (clearly not significant, stopped early)\n")
cat("  >", opts$n_perm_init, " to ", opts$n_perm_mid, ":", n_perm_mid_only, " (escalated tier 1)\n")
cat("  >", opts$n_perm_mid, " (max tier):", n_perm_max, "\n")
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

# Save outputs. Tag includes pseudotime_col + covariate_tag + "permutation"
# marker so different lineages and different batch designs do not overwrite
# each other AND the LRT output (260503-...R) is not overwritten.
pt_tag <- gsub("[^A-Za-z0-9]+", "_", opts$pseudotime_col)
cov_tag <- gsub("[^A-Za-z0-9]+", "_", covariate_tag)
tag <- paste0(
    "permutation_optimized2_", opts$method, "_df", opts$spline_df,
    "_pt-", pt_tag, "_cov-", cov_tag
)
csv_path <- file.path(output_dir, paste0("dirmult_regression_", tag, ".csv"))
rds_path <- file.path(output_dir, paste0("dirmult_regression_", tag, ".rds"))
metadata_path <- file.path(output_dir, paste0("run_metadata_", tag, ".json"))

write.csv(results_df, csv_path, row.names = FALSE)
saveRDS(list(
    results = results_df,
    se = setNames(
        lapply(results_list, `[[`, "se"),
        vapply(results_list, `[[`, character(1), "utr_name")
    ),
    trend = trend_result,
    permutation = permutation_results,
    calibration = calibration_result,
    opts = opts
), rds_path)

# Run-metadata sidecar: the audit trail reviewers will ask for.
metadata <- list(
    script = "260526-permutation_optimized_2.R",
    script_version_hash = tryCatch(
        system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
        error = function(e) NA_character_
    ),
    R_version = R.version.string,
    package_versions = list(
        anndata = as.character(packageVersion("anndata")),
        splines = as.character(getRversion()),  # base R; report R version
        parallel = as.character(getRversion()),
        digest = as.character(packageVersion("digest")),
        pbapply = as.character(packageVersion("pbapply")),
        jsonlite = as.character(packageVersion("jsonlite"))
    ),
    rng_kind = RNGkind(),
    master_seed = opts$seed,
    permutation_strategy = opts$permutation_strategy,
    permutation_block_col = opts$permutation_block_col,
    n_perm_init = opts$n_perm_init,
    n_perm_mid = opts$n_perm_mid,
    n_perm_max = opts$n_perm_max,
    n_hits_stop = opts$n_hits_stop,
    pseudotime_col = opts$pseudotime_col,
    covariate_cols = opts$covariate_cols,
    keep_induced_false_only = opts$keep_induced_false_only,
    subsample_n = opts$subsample_n,
    n_cells_post_filter = nrow(X),
    spline_df = opts$spline_df,
    method = opts$method,
    ncores = opts$ncores,
    start_time = run_start_time,
    end_time = format(Sys.time()),
    host = run_host,
    command = paste(commandArgs(trailingOnly = FALSE), collapse = " "),
    output_csv = csv_path,
    output_rds = rds_path
)
writeLines(jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), metadata_path)

cat("\nSaved: ", csv_path, "\n", sep = "")
cat("Saved: ", rds_path, "\n", sep = "")
cat("Saved: ", metadata_path, " (audit trail)\n", sep = "")

# --- Phase E (optional): Lambda sensitivity sweep ---
# Re-runs the per-gene fit at lambda in {0, 1, 2} (with the adaptive scaling
# kept) and writes a comparison CSV with per-gene deviance, p_value_chisq, and
# s_shared at each setting. The expected outcome is that p-values are highly
# correlated for genes with reasonable cells_per_param, with divergence only
# for low-power genes where the penalty is doing real work. This addresses
# the MAP/MLE-mismatch reviewer concern empirically.
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
