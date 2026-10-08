#!/usr/bin/env Rscript

# =============================================================================
# DRIMSeq treatment-effect LRT on (sample x pseudotime-bin) pseudobulk.
# PRIMARY SCREEN: MAIN EFFECT ONLY (no induced:pseudotime interaction term).
#
# Companion script:
#   drimseq_treatment_interaction_pseudobulk.R  -- secondary screen with the
#   interaction term included; use when you want to discover hits whose
#   treatment effect varies along pseudotime but averages out.
#
# Knoblich brain-organoid 3'UTR pipeline, treatment arm.
#
# Per gene per lineage, one LRT is computed:
#
#   full      ~ bs(pt, df) + super_batch + induced
#   null      ~ bs(pt, df) + super_batch
#
# The induced contrast is tested holding pseudotime fixed (via the spline).
# This is the same model the previous 3-test version called the "main" test.
# Why drop the interaction? With n=164 pseudobulk units, the full design that
# also includes induced:bs(pt) has 14 columns -> 43 free parameters per gene
# at K=4 (3.8 obs/param, below the 5x rule of thumb). dmFit's BFGS optimiser
# silently fails to converge for K>=4, dropping ~33% of post-filter genes.
# Removing the interaction reduces to 9 columns -> 28 free parameters at K=4
# (5.9 obs/param), and gene-wise precision becomes reliable. See
# 260511-drimseq_treatment_methods.md Section 7-11 for the full rationale.
#
# Per-gene precision strategy: gene-wise with empirical-Bayes shrinkage as
# DRIMSeq's default. A small residual fraction (typically K>=7) still gets
# non-finite gene-wise precision; those genes have their precision overwritten
# with the common precision. This is a much lighter-touch correction than the
# force-everything-to-common used in the interaction script.
#
# A cluster bootstrap by sample_id is available as a supplementary sensitivity
# analysis (off by default; n=22 samples is borderline).
#
# Companion methods document:
#   260511-drimseq_treatment_methods.md
#
# Usage (one lineage per invocation; the shell wrapper loops):
#
#   Rscript drimseq_treatment_main_pseudobulk.R \
#       --input  /path/to/merged_*.h5ad \
#       --lineage inhibitory \
#       --output_dir /path/to/results/inhibitory/main_only \
#       --spline_df 5 --n_bins 8 --min_cells_per_unit 10 \
#       --bootstrap_b 0 --ncores 16
#
# Test mode:
#   Rscript drimseq_treatment_main_pseudobulk.R --test
# =============================================================================

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
    library(rlang)
    library(Matrix)
    library(splines)
})

reticulate::use_condaenv("R_downstream", required = TRUE)
box::use(
    DRIMSeq[
        dmDSdata,
        dmFilter,
        dmPrecision,
        dmFit,
        dmTest,
        results,
        samples,
        counts,
        common_precision
    ],
    dplyr[
        arrange,
        mutate,
        select,
        filter,
        `%>%`,
        bind_rows,
        left_join,
        group_by,
        summarise,
        n_distinct,
        n
    ],
    tibble[tibble, as_tibble]
)

# =============================================================================
# Section 0: Argument parsing
# =============================================================================

parse_args <- function(argv) {
    defaults <- list(
        input = paste0(DATA_ROOT, "/3utr/organoids/knoblich/utrcalling/merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad"),
        lineage = "inhibitory",
        output_dir = NULL,
        require_targeting_guide = FALSE,
        control_guide = "Control2",
        balanced_batches_only = FALSE,
        spline_df = 5L,
        n_bins = 8L,
        min_cells_per_unit = 10L,
        min_samps_gene_expr = 5L,
        min_gene_expr = 50L,
        min_samps_feature_expr = 5L,
        min_feature_expr = 20L,
        min_units_per_cond = 3L,
        min_cond_expr = 10L,
        bootstrap_b = 0L,
        bootstrap_seed = 2026L,
        permute_labels = 0L,
        permute_seed = 2027L,
        seed = 42L,
        ncores = max(1L, parallel::detectCores() - 1L),
        stability_top_n = 100L,
        stability_top_2n = 200L,
        test = FALSE
    )
    i <- 1L
    while (i <= length(argv)) {
        a <- argv[i]
        if (a == "--test") {
            defaults$test <- TRUE
            i <- i + 1L
            next
        }
        if (startsWith(a, "--") && i < length(argv)) {
            key <- sub("^--", "", a)
            val <- argv[i + 1L]
            if (!(key %in% names(defaults))) {
                stop("Unknown argument: ", a)
            }
            cur <- defaults[[key]]
            defaults[[key]] <- if (!is.null(cur) && is.integer(cur)) {
                as.integer(val)
            } else {
                val
            }
            i <- i + 2L
        } else {
            stop("Bad argument: ", a)
        }
    }
    if (!defaults$test && is.null(defaults$output_dir)) {
        stop("--output_dir is required when not in --test mode")
    }
    defaults
}

# =============================================================================
# Section 1: Helper functions (shared with interaction script -- keep in sync)
# =============================================================================

pseudobulk_by_group <- function(X_cells_by_features, groups) {
    stopifnot(nrow(X_cells_by_features) == length(groups))
    group_factor <- droplevels(as.factor(groups))
    G <- Matrix::fac2sparse(group_factor)
    pb <- G %*% X_cells_by_features
    pb <- Matrix::t(pb)
    colnames(pb) <- levels(group_factor)
    pb
}

derive_super_batch <- function(sample_id_vec) {
    substr(as.character(sample_id_vec), 1L, 3L)
}

# Permute a label vector WITHIN strata (preserves each stratum's label composition).
# Used by the sample-level treatment-label permutation null: shuffling `induced` across
# samples within each super_batch keeps the super_batch x induced design structure fixed
# (super_batch 137, being all-induced, is invariant) while destroying the treatment signal.
permute_within_strata <- function(labels, strata) {
    out <- labels
    for (s in unique(strata)) {
        idx <- which(strata == s)
        if (length(idx) > 1L) out[idx] <- sample(labels[idx])
    }
    out
}

to_dgCMatrix_double <- function(X) {
    if (inherits(X, "RsparseMatrix")) {
        tmp <- as(X, "TsparseMatrix")
        return(Matrix::sparseMatrix(
            i = tmp@i + 1L,
            j = tmp@j + 1L,
            x = as.numeric(tmp@x),
            dims = dim(tmp),
            dimnames = dimnames(tmp)
        ))
    }
    if (!inherits(X, "dgCMatrix")) {
        X <- as(X, "CsparseMatrix")
    }
    if (typeof(X@x) != "double") {
        X@x <- as.numeric(X@x)
    }
    X
}

condition_balance_prefilter <- function(
    counts_df,
    samples_df,
    min_units_per_cond,
    min_cond_expr
) {
    sample_cols <- samples_df$sample_id
    stopifnot(all(sample_cols %in% colnames(counts_df)))
    is_true_chr <- as.character(samples_df$induced) %in%
        c("TRUE", "True", "true", "T", "1")
    is_false_chr <- as.character(samples_df$induced) %in%
        c("FALSE", "False", "false", "F", "0")
    units_T_cols <- sample_cols[is_true_chr]
    units_F_cols <- sample_cols[is_false_chr]
    if (length(units_T_cols) == 0L || length(units_F_cols) == 0L) {
        stop("condition_balance_prefilter: one of the induced groups is empty.")
    }
    ct_T <- as.matrix(counts_df[, units_T_cols, drop = FALSE])
    ct_F <- as.matrix(counts_df[, units_F_cols, drop = FALSE])
    storage.mode(ct_T) <- "double"
    storage.mode(ct_F) <- "double"
    gid <- counts_df$gene_id
    gene_T <- rowsum(ct_T, group = gid, reorder = FALSE)
    gene_F <- rowsum(ct_F, group = gid, reorder = FALSE)
    n_units_T <- rowSums(gene_T >= min_cond_expr)
    n_units_F <- rowSums(gene_F >= min_cond_expr)
    keep_genes <- (n_units_T >= min_units_per_cond) &
        (n_units_F >= min_units_per_cond)
    kept_gene_ids <- rownames(gene_T)[keep_genes]
    keep_rows <- counts_df$gene_id %in% kept_gene_ids
    list(
        counts_df_kept = counts_df[keep_rows, , drop = FALSE],
        n_genes_pre = nrow(gene_T),
        n_genes_post = sum(keep_genes),
        n_features_pre = nrow(counts_df),
        n_features_post = sum(keep_rows)
    )
}

quantile_bin_pseudotime <- function(pt, n_bins) {
    qs <- stats::quantile(
        pt,
        probs = seq(0, 1, length.out = n_bins + 1L),
        na.rm = TRUE,
        names = FALSE
    )
    qs <- unique(qs)
    if (length(qs) < 2L) {
        stop("Pseudotime has fewer than 2 unique values; cannot bin.")
    }
    bin <- cut(pt, breaks = qs, include.lowest = TRUE, labels = FALSE)
    bin_mid <- vapply(
        seq_len(max(bin, na.rm = TRUE)),
        function(b) {
            mean(pt[bin == b], na.rm = TRUE)
        },
        numeric(1)
    )
    list(bin = bin, bin_mid = bin_mid)
}

# Simpler than the interaction script: full has only the additive `induced`
# coefficient (no induced:bs_pt term). LRT df per compositional dimension = 1.
build_design_matrices <- function(samples_table, spline_df) {
    s <- samples_table
    s$super_batch <- factor(as.character(s$super_batch))
    s$induced <- factor(as.character(s$induced), levels = c("FALSE", "TRUE"))
    bs_basis <- splines::bs(s$pt_mid, df = spline_df)
    s$bs_pt <- bs_basis
    list(
        full = model.matrix(~ bs_pt + super_batch + induced, data = s),
        null = model.matrix(~ bs_pt + super_batch, data = s)
    )
}

# Delta_psi computation (verbatim from the interaction script -- shared
# definition of the treatment effect-size signal).
compute_delta_psi <- function(counts_d, samples_df) {
    pb_cols_raw <- samples_df$sample_id
    if (all(pb_cols_raw %in% colnames(counts_d))) {
        pb_cols <- pb_cols_raw
    } else {
        pb_cols <- make.names(pb_cols_raw)
        stopifnot(all(pb_cols %in% colnames(counts_d)))
    }
    feature_meta <- counts_d[, c("gene_id", "feature_id")]
    counts_mat <- as.matrix(counts_d[, pb_cols, drop = FALSE])
    storage.mode(counts_mat) <- "double"
    gene_id <- feature_meta$gene_id
    gene_unique <- unique(gene_id)
    psi_T_vec <- numeric(length(gene_unique))
    psi_F_vec <- numeric(length(gene_unique))
    n_units_T_vec <- integer(length(gene_unique))
    n_units_F_vec <- integer(length(gene_unique))
    n_cells_T_vec <- numeric(length(gene_unique))
    n_cells_F_vec <- numeric(length(gene_unique))
    feat_pi_T <- numeric(nrow(feature_meta))
    feat_pi_F <- numeric(nrow(feature_meta))
    induced_chr <- as.character(samples_df$induced)
    is_T <- induced_chr %in% c("TRUE", "True", "true", "T", "1")
    is_F <- !is_T
    n_cells <- as.numeric(samples_df$n_cells)
    w_T <- n_cells * is_T
    w_F <- n_cells * is_F
    for (gi in seq_along(gene_unique)) {
        g <- gene_unique[gi]
        rows <- which(gene_id == g)
        K <- length(rows)
        if (K < 2L) {
            psi_T_vec[gi] <- NA_real_
            psi_F_vec[gi] <- NA_real_
            n_units_T_vec[gi] <- 0L
            n_units_F_vec[gi] <- 0L
            n_cells_T_vec[gi] <- NA_real_
            n_cells_F_vec[gi] <- NA_real_
            feat_pi_T[rows] <- NA_real_
            feat_pi_F[rows] <- NA_real_
            next
        }
        rank_norm <- (seq_len(K) - 1) / (K - 1)
        Y_g <- counts_mat[rows, , drop = FALSE]
        unit_totals <- colSums(Y_g)
        keep_u <- unit_totals > 0
        if (sum(keep_u) == 0L) {
            psi_T_vec[gi] <- NA_real_
            psi_F_vec[gi] <- NA_real_
            n_units_T_vec[gi] <- 0L
            n_units_F_vec[gi] <- 0L
            n_cells_T_vec[gi] <- NA_real_
            n_cells_F_vec[gi] <- NA_real_
            feat_pi_T[rows] <- NA_real_
            feat_pi_F[rows] <- NA_real_
            next
        }
        Y_keep <- Y_g[, keep_u, drop = FALSE]
        tot_keep <- unit_totals[keep_u]
        prop_uk <- sweep(Y_keep, 2L, tot_keep, "/")
        psi_u <- as.numeric(crossprod(rank_norm, prop_uk))
        is_T_k <- is_T[keep_u]
        is_F_k <- is_F[keep_u]
        w_T_k <- w_T[keep_u]
        w_F_k <- w_F[keep_u]
        sum_w_T_k <- sum(w_T_k)
        sum_w_F_k <- sum(w_F_k)
        if (sum_w_T_k > 0) {
            psi_T_vec[gi] <- sum(psi_u * w_T_k) / sum_w_T_k
            feat_pi_T[rows] <- as.numeric(prop_uk %*% (w_T_k / sum_w_T_k))
        } else {
            psi_T_vec[gi] <- NA_real_
            feat_pi_T[rows] <- NA_real_
        }
        if (sum_w_F_k > 0) {
            psi_F_vec[gi] <- sum(psi_u * w_F_k) / sum_w_F_k
            feat_pi_F[rows] <- as.numeric(prop_uk %*% (w_F_k / sum_w_F_k))
        } else {
            psi_F_vec[gi] <- NA_real_
            feat_pi_F[rows] <- NA_real_
        }
        n_units_T_vec[gi] <- sum(is_T_k)
        n_units_F_vec[gi] <- sum(is_F_k)
        n_cells_T_vec[gi] <- sum(w_T_k)
        n_cells_F_vec[gi] <- sum(w_F_k)
    }
    delta_psi_gene <- tibble::tibble(
        gene_id = gene_unique,
        psi_induced_True = psi_T_vec,
        psi_induced_False = psi_F_vec,
        delta_psi_treatment = psi_T_vec - psi_F_vec,
        n_units_True = n_units_T_vec,
        n_units_False = n_units_F_vec,
        n_cells_True = n_cells_T_vec,
        n_cells_False = n_cells_F_vec
    )
    delta_psi_feature <- tibble::tibble(
        gene_id = feature_meta$gene_id,
        feature_id = feature_meta$feature_id,
        pi_induced_True = feat_pi_T,
        pi_induced_False = feat_pi_F,
        delta_pi_feature = feat_pi_T - feat_pi_F
    )
    list(gene = delta_psi_gene, feature = delta_psi_feature)
}

# =============================================================================
# Section 2: DRIMSeq pipeline wrapper -- single LRT, gene-wise precision with
# NaN repair.
# =============================================================================
# Precision strategy: leave gene-wise precision as DRIMSeq estimates it, but
# replace NaN/Inf/<=0 entries (typically K>=7 in this dataset) with the
# common precision. This is the lighter-touch repair appropriate for the
# simpler "main only" design, where gene-wise estimates are well-behaved for
# the bulk of the data (verified 2026-05-13 on full inhibitory: 5,291/5,441
# valid p-values, 76 hits at adj_p<0.01).
run_drimseq_pipeline <- function(
    d,
    design_full,
    design_null,
    bpp,
    add_uniform = TRUE,
    verbose = TRUE
) {
    if (verbose) {
        message("  dmPrecision (genewise)...")
    }
    used_genewise <- TRUE
    d <- tryCatch(
        DRIMSeq::dmPrecision(
            d,
            design = design_full,
            add_uniform = add_uniform,
            BPPARAM = bpp
        ),
        error = function(e) {
            if (verbose) {
                message("  genewise precision failed: ", conditionMessage(e))
                message("  retrying with common precision only...")
            }
            used_genewise <<- FALSE
            DRIMSeq::dmPrecision(
                d,
                design = design_full,
                genewise_precision = FALSE,
                add_uniform = add_uniform,
                BPPARAM = bpp
            )
        }
    )

    # NaN repair: replace non-finite / zero / negative gene-wise precisions
    # with the common precision. Without this, dmFit silently produces
    # non-finite log-likelihoods for the affected genes and dmTest reports
    # NA p-values. The slot is a named numeric vector; we manipulate it
    # directly because the public setter doesn't roundtrip names cleanly.
    n_repaired <- 0L
    repair_outcome <- tryCatch(
        {
            if (!used_genewise) {
                if (verbose) {
                    message(
                        "  precision: dmPrecision already ran with genewise=FALSE; skipping NaN repair"
                    )
                }
                "skipped_common"
            } else {
                gw_vec <- d@genewise_precision
                cp <- as.numeric(d@common_precision)[1]
                stopifnot(
                    is.numeric(gw_vec),
                    length(gw_vec) > 0L,
                    is.finite(cp),
                    cp > 0
                )
                bad <- !is.finite(gw_vec) | gw_vec <= 0
                n_repaired <<- as.integer(sum(bad))
                if (n_repaired > 0L) {
                    gw_new <- gw_vec
                    gw_new[bad] <- cp
                    d@genewise_precision <- gw_new
                    if (verbose) {
                        message(
                            "  precision NaN repair: ",
                            n_repaired,
                            "/",
                            length(gw_vec),
                            " gene-wise precisions had NaN/Inf/<=0; ",
                            "replaced with common_precision (",
                            round(cp, 4),
                            ")"
                        )
                    }
                    # Sanity check
                    gw_check_df <- DRIMSeq::genewise_precision(d)
                    stopifnot(
                        is.data.frame(gw_check_df),
                        nrow(gw_check_df) == length(gw_new),
                        all(is.finite(gw_check_df$genewise_precision)),
                        all(gw_check_df$genewise_precision > 0)
                    )
                    "repaired"
                } else {
                    if (verbose) {
                        message(
                            "  precision: 0 genes needed NaN repair (all ",
                            length(gw_vec),
                            " have finite positive precision)"
                        )
                    }
                    "ok"
                }
            }
        },
        error = function(e) {
            if (verbose) {
                message(
                    "  WARNING: precision NaN repair block failed (",
                    conditionMessage(e),
                    "); proceeding with whatever dmPrecision produced"
                )
            }
            "errored"
        }
    )

    if (verbose) {
        message("  dmFit...")
    }
    d <- DRIMSeq::dmFit(d, design = design_full, BPPARAM = bpp)

    if (verbose) {
        message("  dmTest (treatment main)...")
    }
    d_test <- DRIMSeq::dmTest(d, design = design_null, BPPARAM = bpp)

    list(
        d_fit = d,
        gene = DRIMSeq::results(d_test, level = "gene") |>
            tibble::as_tibble() |>
            arrange(pvalue),
        feature = DRIMSeq::results(d_test, level = "feature") |>
            tibble::as_tibble() |>
            arrange(pvalue),
        common_precision = tryCatch(
            as.numeric(DRIMSeq::common_precision(d)),
            error = function(e) NA_real_
        ),
        n_genes_precision_repaired = n_repaired
    )
}

# =============================================================================
# Section 3: Cluster bootstrap by sample_id (single-test version)
# =============================================================================
cluster_bootstrap <- function(
    counts_df,
    pb_samples_drim,
    spline_df,
    dm_filter_args,
    bpp,
    B,
    seed,
    verbose = TRUE
) {
    stopifnot(B > 0L)
    raw_ids <- pb_samples_drim$sample_id
    mangled_ids <- make.names(raw_ids)
    cn <- colnames(counts_df)
    if (!all(raw_ids %in% cn) && all(mangled_ids %in% cn)) {
        idx <- match(mangled_ids, cn)
        colnames(counts_df)[idx] <- raw_ids
    }
    unique_samples <- unique(pb_samples_drim$sample_id_root)
    n_samples <- length(unique_samples)
    if (verbose) {
        message("Cluster bootstrap: ", n_samples, " unique samples, B = ", B)
        if (n_samples < 30L) {
            message(
                "  WARN: n_samples (",
                n_samples,
                ") below the standard cluster-bootstrap reliability threshold of 30."
            )
        }
    }
    set.seed(seed)
    boot_seeds <- sample.int(.Machine$integer.max, B)
    gene_id_observed <- unique(counts_df$gene_id)
    n_genes <- length(gene_id_observed)
    pvalue_mat <- matrix(
        NA_real_,
        nrow = n_genes,
        ncol = B,
        dimnames = list(gene_id_observed, NULL)
    )
    for (b in seq_len(B)) {
        if (verbose) {
            cat(sprintf("  bootstrap %d / %d ... ", b, B))
        }
        set.seed(boot_seeds[b])
        chosen <- sample(unique_samples, size = n_samples, replace = TRUE)
        units_per_sample <- split(
            pb_samples_drim,
            pb_samples_drim$sample_id_root
        )
        boot_count_cols <- vector("list", length(chosen))
        boot_samples_list <- vector("list", length(chosen))
        for (j in seq_along(chosen)) {
            samp <- chosen[j]
            unit_meta <- units_per_sample[[samp]]
            unit_cols <- unit_meta$sample_id
            sub_cnt <- counts_df[, unit_cols, drop = FALSE]
            new_names <- paste0(unit_cols, "__rep", j)
            colnames(sub_cnt) <- new_names
            unit_meta$sample_id <- new_names
            boot_count_cols[[j]] <- sub_cnt
            boot_samples_list[[j]] <- unit_meta
        }
        boot_counts <- do.call(cbind, boot_count_cols)
        boot_counts_df <- data.frame(
            gene_id = counts_df$gene_id,
            feature_id = counts_df$feature_id,
            boot_counts,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        boot_samples_df <- do.call(rbind, boot_samples_list)
        rownames(boot_samples_df) <- NULL
        stopifnot(all(colnames(boot_counts) == boot_samples_df$sample_id))
        d_b <- tryCatch(
            DRIMSeq::dmDSdata(
                counts = boot_counts_df,
                samples = as.data.frame(boot_samples_df)
            ),
            error = function(e) NULL
        )
        if (is.null(d_b)) {
            if (verbose) {
                cat("dmDSdata failed\n")
            }
            next
        }
        d_b <- tryCatch(
            do.call(DRIMSeq::dmFilter, c(list(d_b), dm_filter_args)),
            error = function(e) NULL
        )
        if (is.null(d_b) || nrow(DRIMSeq::counts(d_b)) == 0L) {
            if (verbose) {
                cat("filter dropped all\n")
            }
            next
        }
        designs <- tryCatch(
            build_design_matrices(DRIMSeq::samples(d_b), spline_df),
            error = function(e) NULL
        )
        if (is.null(designs)) {
            if (verbose) {
                cat("design failed\n")
            }
            next
        }
        if (qr(designs$full)$rank < ncol(designs$full)) {
            if (verbose) {
                cat("rank-deficient\n")
            }
            next
        }
        out <- tryCatch(
            run_drimseq_pipeline(
                d_b,
                designs$full,
                designs$null,
                bpp = bpp,
                verbose = FALSE
            ),
            error = function(e) {
                if (verbose) {
                    cat("pipeline failed: ", conditionMessage(e), "\n")
                }
                NULL
            }
        )
        if (is.null(out)) {
            next
        }
        idx <- match(out$gene$gene_id, gene_id_observed)
        keep_idx <- !is.na(idx)
        pvalue_mat[idx[keep_idx], b] <- out$gene$pvalue[keep_idx]
        if (verbose) cat("done\n")
    }
    pvalue_mat
}

# =============================================================================
# Section 3b: Sample-level treatment-label permutation null
# =============================================================================
# Assumption-light calibration of the pseudobulk DRIMSeq LRT for the `induced`
# contrast, and the proper analogue of the within-sample pseudotime permutation
# used for Figure 1 / the UTR pseudotime screen. Because `induced` is CONSTANT
# within sample, the only valid shuffle is at the SAMPLE level; we permute it
# WITHIN super_batch (stratified) so the super_batch nuisance structure is
# preserved (super_batch 137 is all-induced and therefore invariant) and only
# the treatment signal is destroyed. The post-filter gene set (counts_post) is
# held fixed; we relabel the samples table and re-run dmPrecision/dmFit/dmTest.
# Returns genes x B matrices of permutation-null gene-level p-values and LR stats.
#
# Note: with ~14 effective samples and the super_batch constraint, the stratified
# permutation space is small (~hundreds of arrangements), so this CALIBRATES the
# chi-squared DRIMSeq p-values (per-gene permutation-p floor 1/(B+1); aggregate
# hit-count null) rather than replacing them - exactly the role permutation plays
# for the Figure 1 LRT. See TODO.md.
permute_labels_null <- function(
    counts_post,
    pb_samples_drim,
    spline_df,
    dm_filter_args,
    bpp,
    B,
    seed,
    verbose = TRUE
) {
    stopifnot(B > 0L)
    su <- unique(pb_samples_drim[, c(
        "sample_id_root",
        "induced",
        "super_batch"
    )])
    rownames(su) <- su$sample_id_root
    gene_ids <- unique(counts_post$gene_id)
    ng <- length(gene_ids)
    pmat <- matrix(NA_real_, ng, B, dimnames = list(gene_ids, NULL))
    lrmat <- matrix(NA_real_, ng, B, dimnames = list(gene_ids, NULL))
    n_perm_sb <- sum(vapply(
        unique(su$super_batch),
        function(s) length(unique(su$induced[su$super_batch == s])) > 1L,
        logical(1)
    ))
    if (verbose) {
        message(
            "Label permutation: ",
            nrow(su),
            " samples; ",
            n_perm_sb,
            " super_batch(es) carry both labels; B = ",
            B
        )
    }
    set.seed(seed)
    perm_seeds <- sample.int(.Machine$integer.max, B)
    for (b in seq_len(B)) {
        if (verbose) {
            cat(sprintf("  permutation %d / %d ... ", b, B))
        }
        set.seed(perm_seeds[b])
        su_perm_lab <- permute_within_strata(su$induced, su$super_batch)
        names(su_perm_lab) <- su$sample_id_root
        samples_perm <- pb_samples_drim
        samples_perm$induced <- su_perm_lab[samples_perm$sample_id_root]
        d_b <- tryCatch(
            DRIMSeq::dmDSdata(
                counts = counts_post,
                samples = as.data.frame(samples_perm)
            ),
            error = function(e) NULL
        )
        if (is.null(d_b)) {
            if (verbose) {
                cat("dmDSdata failed\n")
            }
            next
        }
        d_b <- tryCatch(
            do.call(DRIMSeq::dmFilter, c(list(d_b), dm_filter_args)),
            error = function(e) NULL
        )
        if (is.null(d_b) || nrow(DRIMSeq::counts(d_b)) == 0L) {
            if (verbose) {
                cat("filter dropped all\n")
            }
            next
        }
        designs <- tryCatch(
            build_design_matrices(DRIMSeq::samples(d_b), spline_df),
            error = function(e) NULL
        )
        if (is.null(designs) || qr(designs$full)$rank < ncol(designs$full)) {
            if (verbose) {
                cat("rank-deficient\n")
            }
            next
        }
        out <- tryCatch(
            run_drimseq_pipeline(
                d_b,
                designs$full,
                designs$null,
                bpp = bpp,
                verbose = FALSE
            ),
            error = function(e) {
                if (verbose) {
                    cat("pipeline failed\n")
                }
                NULL
            }
        )
        if (is.null(out)) {
            next
        }
        idx <- match(out$gene$gene_id, gene_ids)
        keep <- !is.na(idx)
        pmat[idx[keep], b] <- out$gene$pvalue[keep]
        lrmat[idx[keep], b] <- out$gene$lr[keep]
        if (verbose) cat("done\n")
    }
    list(pmat = pmat, lrmat = lrmat)
}

# =============================================================================
# Section 4: Test mode (--test)
# =============================================================================

run_tests <- function() {
    cat("=== drimseq_treatment_main_pseudobulk: verification tests ===\n\n")
    all_pass <- TRUE

    cat("Test 1: pseudobulk_by_group correctness... ")
    set.seed(42)
    n_cells <- 100L
    n_feat <- 7L
    dense_in <- matrix(as.numeric(rpois(n_cells * n_feat, 3)), nrow = n_cells)
    cells_x_feat <- as(Matrix::Matrix(dense_in, sparse = TRUE), "CsparseMatrix")
    grp <- sample(letters[1:4], n_cells, replace = TRUE)
    pb <- pseudobulk_by_group(cells_x_feat, grp)
    expected <- vapply(
        sort(unique(grp)),
        function(g) {
            colSums(as.matrix(cells_x_feat[grp == g, , drop = FALSE]))
        },
        numeric(n_feat)
    )
    rownames(expected) <- NULL
    if (max(abs(as.matrix(pb) - expected)) < 1e-10) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    cat("Test 2: derive_super_batch... ")
    sb <- derive_super_batch(c("137001", "206alpha", "213xx", "9999"))
    if (identical(sb, c("137", "206", "213", "999"))) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    cat("Test 3: quantile_bin_pseudotime... ")
    bres <- quantile_bin_pseudotime(runif(1000), 8L)
    tab <- table(bres$bin)
    if (length(tab) == 8L && max(tab) - min(tab) <= 2L) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    cat("Test 4: design-matrix df arithmetic (main-only)... ")
    set.seed(42)
    sdf <- 5L
    n_units <- 100L
    samples_table <- data.frame(
        sample_id = sprintf("u%03d", seq_len(n_units)),
        pt_mid = sort(runif(n_units)),
        super_batch = sample(c("137", "206", "213"), n_units, replace = TRUE),
        induced = sample(c("FALSE", "TRUE"), n_units, replace = TRUE)
    )
    designs <- build_design_matrices(samples_table, spline_df = sdf)
    # full: 1 (int) + sdf + 2 (sb-1) + 1 (induced) = 4 + sdf
    # null: 1 + sdf + 2 = 3 + sdf
    # df per K-1 = 1
    df_main <- ncol(designs$full) - ncol(designs$null)
    if (
        ncol(designs$full) == 4L + sdf &&
            ncol(designs$null) == 3L + sdf &&
            df_main == 1L
    ) {
        cat(
            "PASS (full=",
            ncol(designs$full),
            " null=",
            ncol(designs$null),
            " df=",
            df_main,
            ")\n",
            sep = ""
        )
    } else {
        cat(
            "FAIL (full=",
            ncol(designs$full),
            " null=",
            ncol(designs$null),
            " df=",
            df_main,
            ")\n",
            sep = ""
        )
        all_pass <- FALSE
    }

    cat("Test 5: compute_delta_psi correctness... ")
    counts_d <- data.frame(
        gene_id = c("A", "A", "A", "B", "B"),
        feature_id = c("A1", "A2", "A3", "B1", "B2"),
        u1 = c(80, 10, 10, 50, 50),
        u2 = c(10, 10, 80, 40, 60),
        u3 = c(10, 20, 70, 30, 70),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    samples_df <- data.frame(
        sample_id = c("u1", "u2", "u3"),
        induced = c("FALSE", "TRUE", "TRUE"),
        n_cells = c(10, 20, 30),
        stringsAsFactors = FALSE
    )
    dp <- compute_delta_psi(counts_d, samples_df)
    a_row <- dp$gene[dp$gene$gene_id == "A", ]
    if (abs(a_row$delta_psi_treatment - 0.67) < 1e-8) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    cat("Test 6: induced reference level... ")
    s2 <- data.frame(
        sample_id = sprintf("u%02d", seq_len(20L)),
        pt_mid = seq(0.05, 0.95, length.out = 20L),
        super_batch = rep(c("206", "213"), each = 10L),
        induced = rep(c("FALSE", "TRUE"), times = 10L),
        stringsAsFactors = FALSE
    )
    designs2 <- build_design_matrices(s2, spline_df = 3L)
    if ("inducedTRUE" %in% colnames(designs2$full)) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    cat("Test 7: end-to-end pipeline runs under H0 (no interaction)... ")
    set.seed(42L)
    n_genes <- 200L
    n_bins_t <- 8L
    pi_true <- c(0.3, 0.3, 0.4)
    K_t <- length(pi_true)
    s_true <- 50
    sample_meta <- data.frame(
        sample_id_root = sprintf("S%02d", seq_len(22L)),
        super_batch = c(rep("137", 8L), rep("213", 4L), rep("206", 10L)),
        induced = c(
            rep("TRUE", 8L),
            rep("TRUE", 3L),
            rep("FALSE", 1L),
            rep("TRUE", 3L),
            rep("FALSE", 7L)
        ),
        stringsAsFactors = FALSE
    )
    units <- expand.grid(
        sample_id_root = sample_meta$sample_id_root,
        bin_id = seq_len(n_bins_t),
        stringsAsFactors = FALSE
    )
    units$sample_id <- paste(units$sample_id_root, units$bin_id, sep = "__")
    units <- merge(units, sample_meta, by = "sample_id_root")
    units$pt_mid <- (units$bin_id - 0.5) / n_bins_t
    units$n_cells <- 100L
    n_units_t <- nrow(units)
    counts_list <- vector("list", n_genes)
    alpha <- pi_true * s_true
    for (gi in seq_len(n_genes)) {
        Y <- t(vapply(
            seq_len(n_units_t),
            function(i) {
                p_dir <- rgamma(K_t, alpha, 1)
                p_dir <- p_dir / sum(p_dir)
                as.integer(rmultinom(1, units$n_cells[i], p_dir))
            },
            integer(K_t)
        ))
        df_g <- data.frame(
            gene_id = paste0("G", gi),
            feature_id = paste0("G", gi, "_", seq_len(K_t)),
            t(Y),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        colnames(df_g)[3:(2 + n_units_t)] <- units$sample_id
        counts_list[[gi]] <- df_g
    }
    counts_df <- do.call(rbind, counts_list)
    samples_drim <- data.frame(
        sample_id = units$sample_id,
        sample_id_root = units$sample_id_root,
        bin_id = units$bin_id,
        pt_mid = units$pt_mid,
        induced = units$induced,
        super_batch = units$super_batch,
        n_cells = units$n_cells,
        stringsAsFactors = FALSE
    )
    bpp <- BiocParallel::SerialParam()
    d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples_drim)
    d <- DRIMSeq::dmFilter(
        d,
        min_samps_gene_expr = 5L,
        min_gene_expr = 50L,
        min_samps_feature_expr = 5L,
        min_feature_expr = 20L
    )
    designs <- build_design_matrices(DRIMSeq::samples(d), spline_df = 5L)
    res <- run_drimseq_pipeline(
        d,
        designs$full,
        designs$null,
        bpp = bpp,
        verbose = FALSE
    )
    pvals <- res$gene$pvalue[!is.na(res$gene$pvalue)]
    mean_p <- if (length(pvals) > 0) mean(pvals) else NA_real_
    if (
        length(pvals) / n_genes >= 0.5 &&
            !is.na(mean_p) &&
            mean_p >= 0.25 &&
            mean_p <= 0.75
    ) {
        cat(
            "PASS (n_valid=",
            length(pvals),
            " mean_p=",
            round(mean_p, 3),
            ")\n",
            sep = ""
        )
    } else {
        cat(
            "FAIL (n_valid=",
            length(pvals),
            " mean_p=",
            round(mean_p, 3),
            ")\n",
            sep = ""
        )
        all_pass <- FALSE
    }

    cat(
        "Test 8: permute_within_strata (stratified sample-level label permutation)... "
    )
    lab8 <- c("TRUE", "TRUE", "TRUE", "FALSE", "FALSE", "TRUE", "FALSE")
    str8 <- c("137", "137", "206", "206", "206", "213", "213")
    set.seed(7L)
    p8 <- permute_within_strata(lab8, str8)
    comp_ok <- all(vapply(
        unique(str8),
        function(s) {
            identical(sort(lab8[str8 == s]), sort(p8[str8 == s]))
        },
        logical(1)
    ))
    inv137 <- all(p8[str8 == "137"] == "TRUE") # all-induced stratum is invariant
    if (comp_ok && inv137) {
        cat("PASS\n")
    } else {
        cat("FAIL\n")
        all_pass <- FALSE
    }

    cat(
        "\n=== ",
        if (all_pass) "ALL TESTS PASS" else "TESTS FAILED",
        " ===\n",
        sep = ""
    )
    if (!all_pass) quit(status = 1L)
}

# =============================================================================
# MAIN
# =============================================================================

opts <- parse_args(commandArgs(trailingOnly = TRUE))
if (opts$test) {
    run_tests()
    quit(status = 0L)
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
run_start_time <- format(Sys.time())
bpp <- BiocParallel::MulticoreParam(workers = opts$ncores)
slug <- gsub("[^A-Za-z0-9]+", "_", tolower(opts$lineage))

cat("=== DRIMSeq treatment-effect LRT (MAIN ONLY, no interaction) ===\n")
cat("  input:                ", opts$input, "\n")
cat("  output_dir:           ", opts$output_dir, "\n")
cat("  lineage:              ", opts$lineage, "\n")
cat("  spline_df:            ", opts$spline_df, "\n")
cat("  n_bins:               ", opts$n_bins, "\n")
cat("  min_cells_per_unit:   ", opts$min_cells_per_unit, "\n")
cat("  bootstrap_b:          ", opts$bootstrap_b, "\n")
cat("  ncores:               ", opts$ncores, "\n")
cat("  start:                ", format(Sys.time()), "\n")

# ---- Load + lineage subset ----
cat("\n[Load] Reading h5ad...\n")
adata <- anndata::read_h5ad(opts$input)
required_obs <- c(
    "sample_id",
    "induced",
    paste0("dpt_pseudotime_", opts$lineage)
)
missing_obs <- setdiff(required_obs, colnames(adata$obs))
if (length(missing_obs) > 0L) {
    stop("Missing obs columns: ", paste(missing_obs, collapse = ", "))
}
if (!"gene_name" %in% colnames(adata$var)) {
    stop("Missing var column gene_name")
}

obs_all <- adata$obs
var_all <- adata$var
pt_all <- as.numeric(obs_all[[paste0("dpt_pseudotime_", opts$lineage)]])
sample_id_all <- as.character(obs_all$sample_id)
induced_all <- as.character(obs_all$induced)
super_batch_all <- derive_super_batch(sample_id_all)
keep <- !is.na(pt_all)
# Optional cleanup: restrict the treatment (induced) arm to confirmed targeting guides,
# dropping induced cells that carry the non-targeting control guide or no assigned guide.
# Uninduced cells (the no-perturbation baseline) are kept regardless of guide. OFF by default
# (preserves the published all-induced behaviour); enable with --require_targeting_guide TRUE.
if (isTRUE(as.logical(opts$require_targeting_guide))) {
    if (!"gRNA" %in% colnames(obs_all)) {
        stop("require_targeting_guide=TRUE needs a 'gRNA' obs column")
    }
    gRNA_all <- as.character(obs_all$gRNA)
    is_ind_all <- induced_all %in% c("TRUE", "True", "true", "T", "1")
    drop_ctrl <- is_ind_all & (is.na(gRNA_all) | gRNA_all == opts$control_guide)
    cat(
        "  require_targeting_guide: dropping",
        sum(drop_ctrl & keep),
        "induced cells with control/unassigned guide (treatment = targeting guides only)\n"
    )
    keep <- keep & !drop_ctrl
}
# Optional: restrict to super_batches that carry BOTH induced and uninduced cells
# (drops super_batch 137, which is 100% induced). Removes the induced~super_batch
# confound entirely at the cost of the ~17k induced batch-137 cells. OFF by default.
if (isTRUE(as.logical(opts$balanced_batches_only))) {
    ind_k <- induced_all %in% c("TRUE", "True", "true", "T", "1")
    both_sb <- names(which(tapply(
        ind_k[keep],
        super_batch_all[keep],
        function(z) length(unique(z)) > 1L
    )))
    cat(
        "  balanced_batches_only: keeping super_batches with both conditions:",
        paste(both_sb, collapse = ", "),
        "\n"
    )
    keep <- keep & (super_batch_all %in% both_sb)
}
cat("  cells in lineage:", sum(keep), "/", length(keep), "\n")

X_all <- to_dgCMatrix_double(adata$X)
X_lin <- X_all[keep, , drop = FALSE]
pt_lin <- pt_all[keep]
sample_id_lin <- sample_id_all[keep]
induced_lin <- induced_all[keep]
super_batch_lin <- super_batch_all[keep]
rm(X_all, adata)
gc(verbose = FALSE)

sample_table_lin <- data.frame(
    sample_id = sample_id_lin,
    induced = induced_lin,
    super_batch = super_batch_lin,
    stringsAsFactors = FALSE
)
cat("\nSample composition:\n")
print(table(
    unique(sample_table_lin)$super_batch,
    unique(sample_table_lin)$induced,
    dnn = c("super_batch", "induced")
))

# ---- Pseudobulk ----
cat("\n[Pseudobulk] Quantile-binning pseudotime into", opts$n_bins, "bins...\n")
bres <- quantile_bin_pseudotime(pt_lin, opts$n_bins)
groups <- paste(sample_id_lin, sprintf("bin%02d", bres$bin), sep = "__")
pb_counts <- as.matrix(pseudobulk_by_group(X_lin, groups))
storage.mode(pb_counts) <- "integer"
unit_keys <- colnames(pb_counts)
pb_meta <- tibble(sample_id = unit_keys) |>
    mutate(
        sample_id_root = sub("__bin[0-9]+$", "", sample_id),
        bin_id = as.integer(sub("^.*__bin", "", sample_id))
    )
pb_meta$pt_mid <- bres$bin_mid[pb_meta$bin_id]
sample_lookup <- unique(sample_table_lin)
pb_meta <- left_join(
    pb_meta,
    sample_lookup,
    by = c("sample_id_root" = "sample_id")
)
pb_meta$n_cells <- as.integer(table(groups)[pb_meta$sample_id])
keep_units <- pb_meta$n_cells >= opts$min_cells_per_unit
cat(
    "  units after min_cells_per_unit (>=",
    opts$min_cells_per_unit,
    "): ",
    sum(keep_units),
    "/",
    nrow(pb_meta),
    "\n",
    sep = ""
)
pb_counts <- pb_counts[, keep_units, drop = FALSE]
pb_meta <- pb_meta[keep_units, , drop = FALSE]

manifest_path <- file.path(opts$output_dir, "pseudobulk_units.csv")
write.csv(as.data.frame(pb_meta), manifest_path, row.names = FALSE)
cat("  manifest:", manifest_path, "\n  unit composition:\n")
print(table(
    pb_meta$super_batch,
    pb_meta$induced,
    dnn = c("super_batch", "induced")
))

# ---- DRIMSeq pipeline ----
feature_ids <- rownames(var_all)
stopifnot(!anyDuplicated(feature_ids))
counts_df <- data.frame(
    gene_id = as.character(var_all$gene_name),
    feature_id = feature_ids,
    pb_counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
pb_samples_drim <- as.data.frame(pb_meta)[, c(
    "sample_id",
    "sample_id_root",
    "bin_id",
    "pt_mid",
    "induced",
    "super_batch",
    "n_cells"
)]
stopifnot(all(colnames(pb_counts) == pb_samples_drim$sample_id))

dm_filter_args <- list(
    min_samps_gene_expr = opts$min_samps_gene_expr,
    min_gene_expr = opts$min_gene_expr,
    min_samps_feature_expr = opts$min_samps_feature_expr,
    min_feature_expr = opts$min_feature_expr
)

cat("\n[Pre-filter] Condition-balance pre-filter ...\n")
prefilt <- condition_balance_prefilter(
    counts_df,
    pb_samples_drim,
    opts$min_units_per_cond,
    opts$min_cond_expr
)
cat(
    "  pre :",
    prefilt$n_genes_pre,
    "genes,",
    prefilt$n_features_pre,
    "features\n"
)
cat(
    "  post:",
    prefilt$n_genes_post,
    "genes,",
    prefilt$n_features_post,
    "features\n"
)
if (prefilt$n_features_post == 0L) {
    stop("Pre-filter dropped all genes")
}
counts_df <- prefilt$counts_df_kept

cat("\n[DRIMSeq] dmDSdata + dmFilter...\n")
d <- DRIMSeq::dmDSdata(counts = counts_df, samples = pb_samples_drim)
d <- do.call(DRIMSeq::dmFilter, c(list(d), dm_filter_args))
counts_post <- DRIMSeq::counts(d)
n_genes_pre <- prefilt$n_genes_pre
n_genes_post <- length(unique(counts_post$gene_id))
n_feat_pre <- prefilt$n_features_pre
n_feat_post <- nrow(counts_post)
cat("  after dmFilter:", n_genes_post, "genes,", n_feat_post, "features\n")
if (n_feat_post == 0L) {
    stop("dmFilter dropped all features")
}

cat(
    "\n[DRIMSeq] Building MAIN-ONLY design (spline_df = ",
    opts$spline_df,
    ", no interaction)...\n",
    sep = ""
)
designs <- build_design_matrices(DRIMSeq::samples(d), opts$spline_df)
cat(
    "  full design:",
    ncol(designs$full),
    "cols (rank=",
    qr(designs$full)$rank,
    ")\n",
    sep = ""
)
cat("  null design:", ncol(designs$null), "cols\n")
cat("  LRT df per (K-1):", ncol(designs$full) - ncol(designs$null), "\n")

set.seed(opts$seed)
cat("\n[DRIMSeq] Fit + 1 test (treatment main) on observed data...\n")
obs_out <- run_drimseq_pipeline(
    d,
    designs$full,
    designs$null,
    bpp = bpp,
    verbose = TRUE
)
cat("  common precision =", round(obs_out$common_precision, 4), "\n")
cat(
    "  precision NaN-repaired:",
    obs_out$n_genes_precision_repaired,
    "/",
    n_genes_post,
    "\n",
    sep = ""
)

# ---- Delta_psi ----
cat(
    "\n[Delta_psi] Computing gene-level signed W1 and feature-level Delta_pi...\n"
)
dpsi <- compute_delta_psi(counts_post, DRIMSeq::samples(d))
dpsi_gene_path <- file.path(
    opts$output_dir,
    sprintf("delta_psi_treatment_gene_%s.csv", slug)
)
dpsi_feat_path <- file.path(
    opts$output_dir,
    sprintf("delta_psi_treatment_feature_%s.csv", slug)
)
write.csv(as.data.frame(dpsi$gene), dpsi_gene_path, row.names = FALSE)
write.csv(as.data.frame(dpsi$feature), dpsi_feat_path, row.names = FALSE)

# ---- Output: single gene + feature CSV ----
gene_extras <- counts_post %>%
    group_by(gene_id) %>%
    summarise(n_features_post_filter = n(), .groups = "drop")
gene_extras$n_cells_total <- sum(pb_samples_drim$n_cells)
gene_extras$n_units <- nrow(pb_samples_drim)

gene_df <- obs_out$gene %>%
    left_join(dpsi$gene, by = "gene_id") %>%
    left_join(gene_extras, by = "gene_id") %>%
    mutate(lineage = opts$lineage, test_type = "treatment_main")
feat_df <- obs_out$feature %>%
    left_join(dpsi$feature, by = c("gene_id", "feature_id")) %>%
    mutate(lineage = opts$lineage, test_type = "treatment_main")

gene_csv <- file.path(
    opts$output_dir,
    sprintf("drimseq_gene_%s_treatment_main.csv", slug)
)
feat_csv <- file.path(
    opts$output_dir,
    sprintf("drimseq_feature_diagnostic_%s_treatment_main.csv", slug)
)
write.csv(as.data.frame(gene_df), gene_csv, row.names = FALSE)
write.csv(as.data.frame(feat_df), feat_csv, row.names = FALSE)
cat("\n[Output] ", gene_csv, "\n         ", feat_csv, "\n")

sig_bh05 <- sum(gene_df$adj_pvalue < 0.05, na.rm = TRUE)
sig_bh01 <- sum(gene_df$adj_pvalue < 0.01, na.rm = TRUE)
pvals_test <- gene_df$pvalue[!is.na(gene_df$pvalue)]
ks_p <- if (length(pvals_test) >= 10L) {
    suppressWarnings(stats::ks.test(pvals_test, "punif")$p.value)
} else {
    NA_real_
}

# ---- Cluster bootstrap (optional) ----
if (opts$bootstrap_b > 0L) {
    cat(
        "\n[Bootstrap] Cluster bootstrap by sample_id, B =",
        opts$bootstrap_b,
        "...\n"
    )
    boot_mat <- cluster_bootstrap(
        counts_df = counts_post,
        pb_samples_drim = pb_samples_drim,
        spline_df = opts$spline_df,
        dm_filter_args = dm_filter_args,
        bpp = bpp,
        B = opts$bootstrap_b,
        seed = opts$bootstrap_seed,
        verbose = TRUE
    )

    # Per-gene bootstrap stability for the single test
    gene_ids_boot <- rownames(boot_mat)
    idx <- match(gene_df$gene_id, gene_ids_boot)
    p_obs <- gene_df$pvalue
    n_completed <- rowSums(!is.na(boot_mat))
    n_le_obs <- vapply(
        seq_along(idx),
        function(i) {
            if (is.na(idx[i])) {
                return(NA_integer_)
            }
            bp <- boot_mat[idx[i], ]
            sum(!is.na(bp) & bp <= p_obs[i])
        },
        integer(1)
    )
    pvalue_bootstrap <- (1 + n_le_obs) / (n_completed[idx] + 1)

    # Top-N stability index
    N <- opts$stability_top_n
    twoN <- opts$stability_top_2n
    ordered_obs <- gene_df$gene_id[order(gene_df$pvalue, na.last = TRUE)]
    top_N <- head(ordered_obs, N)
    stab_idx <- vapply(
        seq_len(ncol(boot_mat)),
        function(b) {
            bp <- boot_mat[, b]
            if (all(is.na(bp))) {
                return(NA_real_)
            }
            ordered_b <- names(sort(bp, na.last = TRUE))
            mean(top_N %in% head(ordered_b, twoN))
        },
        numeric(1)
    )
    per_gene_stab <- vapply(
        top_N,
        function(g) {
            if (!g %in% rownames(boot_mat)) {
                return(NA_real_)
            }
            bp <- boot_mat[g, ]
            if (sum(!is.na(bp)) == 0L) {
                return(NA_real_)
            }
            ranks_b <- vapply(
                seq_along(bp),
                function(b) {
                    if (is.na(bp[b])) {
                        return(NA_integer_)
                    }
                    sum(boot_mat[, b] <= bp[b], na.rm = TRUE)
                },
                integer(1)
            )
            mean(ranks_b <= twoN, na.rm = TRUE)
        },
        numeric(1)
    )

    stab_df <- tibble(
        gene_id = gene_df$gene_id,
        pvalue_observed = p_obs,
        adj_pvalue_observed = gene_df$adj_pvalue,
        n_bootstrap_completed = n_completed[idx],
        pvalue_bootstrap = pvalue_bootstrap,
        in_top_N_observed = gene_df$gene_id %in% top_N,
        stability_top_2N = NA_real_
    )
    idx2 <- match(top_N, stab_df$gene_id)
    stab_df$stability_top_2N[idx2] <- per_gene_stab

    stab_path <- file.path(
        opts$output_dir,
        sprintf("bootstrap_stability_gene_%s.csv", slug)
    )
    write.csv(as.data.frame(stab_df), stab_path, row.names = FALSE)
    stab_summary <- tibble(
        bootstrap_rep = seq_along(stab_idx),
        stability_top_N_in_top_2N = stab_idx
    )
    stab_summary_path <- file.path(
        opts$output_dir,
        sprintf("bootstrap_stability_summary_%s.csv", slug)
    )
    write.csv(as.data.frame(stab_summary), stab_summary_path, row.names = FALSE)
    cat(
        "  median top-",
        N,
        "/top-",
        twoN,
        ": ",
        round(median(stab_idx, na.rm = TRUE), 3),
        "\n",
        sep = ""
    )
} else {
    cat("\n[Bootstrap] Skipped (bootstrap_b = 0).\n")
}

# ---- Label-permutation null (optional; assumption-light calibration of the LRT) ----
if (opts$permute_labels > 0L) {
    cat(
        "\n[Permutation] Sample-level treatment-label permutation null, B =",
        opts$permute_labels,
        "...\n"
    )
    perm <- permute_labels_null(
        counts_post = counts_post,
        pb_samples_drim = pb_samples_drim,
        spline_df = opts$spline_df,
        dm_filter_args = dm_filter_args,
        bpp = bpp,
        B = opts$permute_labels,
        seed = opts$permute_seed,
        verbose = TRUE
    )
    obs_lr <- gene_df$lr
    names(obs_lr) <- gene_df$gene_id
    gids <- rownames(perm$lrmat)
    n_done <- rowSums(!is.na(perm$lrmat))
    ge_obs <- vapply(
        seq_along(gids),
        function(i) {
            o <- obs_lr[[gids[i]]]
            v <- perm$lrmat[i, ]
            if (is.null(o) || is.na(o)) {
                return(NA_integer_)
            }
            as.integer(sum(!is.na(v) & v >= o))
        },
        integer(1)
    )
    pvalue_perm <- (1 + ge_obs) / (n_done + 1) # Phipson & Smyth 2010 floor
    perm_gene <- tibble(
        gene_id = gids,
        lr = as.numeric(obs_lr[gids]),
        pvalue = gene_df$pvalue[match(gids, gene_df$gene_id)],
        adj_pvalue = gene_df$adj_pvalue[match(gids, gene_df$gene_id)],
        pvalue_perm = pvalue_perm,
        n_perm_completed = as.integer(n_done)
    )
    perm_path <- file.path(
        opts$output_dir,
        sprintf("permutation_null_gene_%s.csv", slug)
    )
    write.csv(as.data.frame(perm_gene), perm_path, row.names = FALSE)
    hit_null <- colSums(perm$pmat < 0.05, na.rm = TRUE)
    obs_hits_raw <- sum(gene_df$pvalue < 0.05, na.rm = TRUE)
    hc <- tibble(
        permutation = seq_along(hit_null),
        n_genes_p_lt_0p05 = as.integer(hit_null)
    )
    hc_path <- file.path(
        opts$output_dir,
        sprintf("permutation_hitcount_null_%s.csv", slug)
    )
    write.csv(as.data.frame(hc), hc_path, row.names = FALSE)
    cat(sprintf(
        "  observed genes p<0.05: %d ; permutation-null genes p<0.05: mean %.1f, max %d (over %d perms)\n",
        obs_hits_raw,
        mean(hit_null),
        max(hit_null),
        length(hit_null)
    ))
    cat(
        "  [Permutation outputs]",
        perm_path,
        "\n                        ",
        hc_path,
        "\n"
    )
} else {
    cat("\n[Permutation] Skipped (permute_labels = 0).\n")
}

# ---- Summary CSV ----
summary_row <- tibble(
    lineage = opts$lineage,
    n_cells_lineage = sum(keep),
    n_units_pseudobulk = nrow(pb_samples_drim),
    n_genes_pre_filter = n_genes_pre,
    n_genes_post_filter = n_genes_post,
    n_features_pre_filter = n_feat_pre,
    n_features_post_filter = n_feat_post,
    common_precision = obs_out$common_precision,
    n_genes_precision_nan_repaired = obs_out$n_genes_precision_repaired,
    spline_df = opts$spline_df,
    n_bins = opts$n_bins,
    min_cells_per_unit = opts$min_cells_per_unit,
    bootstrap_b = opts$bootstrap_b,
    n_sig_gene_main_bh05 = sig_bh05,
    n_sig_gene_main_bh01 = sig_bh01,
    ks_uniformity_p_main = ks_p,
    start_time = run_start_time,
    end_time = format(Sys.time())
)
summary_path <- file.path(
    opts$output_dir,
    sprintf("drimseq_summary_%s.csv", slug)
)
write.csv(as.data.frame(summary_row), summary_path, row.names = FALSE)
cat("\n[Summary]", summary_path, "\n")
print(as.data.frame(summary_row))
cat("\n=== Done at", format(Sys.time()), "===\n")
