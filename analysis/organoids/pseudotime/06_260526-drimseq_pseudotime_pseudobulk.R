#!/usr/bin/env Rscript

# =============================================================================
# DRIMSeq pseudotime LRT on (sample x pseudotime-bin) pseudobulk.
# Knoblich brain-organoid 3'UTR pipeline, controls-only pseudotime arm.
# Derived from ../perturbation/02_drimseq_treatment_main_pseudobulk.R.
#
# Per gene per lineage (induced=FALSE cells only), one LRT is computed:
#
#   full   ~ batch_number + factor(pt_bin)    (df = 2 for the pt_bin contrast)
#   null   ~ batch_number
#
# batch_number = sample_id = Li et al. 2023 "library", the per-sample
# covariate from their DE model (Y ~ n_UMI + library + cell_type + condition).
# 8 levels in the controls-only subset. The coarser super_batch covariate
# is not used: in the controls-only subset it has only 2 levels, one of
# which is a single sample (see the comment in parse_args()).
#
# pt_bin is a 3-level factor obtained by quantile-binning pseudotime globally
# across all controls cells in the lineage. Each (sample x pt_bin) cell is one
# pseudobulk observation. With 8 controls samples x 3 bins = 24 units per
# lineage, n_eff = 24, enough for a 2-df test of "is there any pt_bin
# effect on UTR usage" after the 7-df adjustment for batch_number. Bins
# below MIN_CELLS_PER_UNIT (default 30) are dropped.
#
# Pseudobulk for pseudotime (same rationale as for the treatment arm):
#   The cell-level LRT in 02_260503-...R assumes per-cell exchangeability under
#   H0. A bootstrap-concordance check on the Knoblich INH lineage showed
#   this assumption fails (concordance 66.7%, Phase D KS p = 0). Pseudobulk
#   by (sample, pt_bin) is the standard comparison without pseudoreplication:
#   each unit is one independent observation, n_eff = n_units. This is the
#   pseudobulk test in the comparison of tests (07_260526-three_way_defensibility).
#
# 3 bins (8 in the treatment script):
#   3 bins is the smallest design that captures monotonic and U-shaped /
#   inverted-U-shaped pseudotime patterns (2 df via factor coding). With 8
#   controls samples this gives 24 units per lineage; pushing to 4 bins
#   would produce a 6-cell unit in INH which is marginal for dispersion
#   estimation.
#
# factor(pt_bin) rather than bs(pt_mid, df=2):
#   Same df, same flexibility, simpler interpretation. The factor reports
#   psi_bin0 / psi_bin1 / psi_bin2 directly; bs() reports spline
#   coefficients which require post-hoc unpacking to be human-readable.
#
# Per-gene precision strategy: gene-wise with empirical-Bayes shrinkage as
# DRIMSeq's default. A small residual fraction (typically K>=7) still gets
# non-finite gene-wise precision; those genes have their precision overwritten
# with the common precision. Same as the treatment-arm script.
#
# Reproducibility:
#   - set.seed(opts$seed) called immediately before each dmPrecision call
#     (dmPrecision uses add_uniform=TRUE which jitters zero-count features).
#   - master_seed recorded in the output CSV + run_metadata.json sidecar.
#   - prec_subset = 1 (forces deterministic precision estimation;
#     dmPrecision is otherwise non-deterministic).
#
# Outputs (gene and feature tables, delta_psi, summary, run_metadata.json)
# are listed in 05_run_pseudobulk_pseudotime.sh.
#
# Usage (one lineage per invocation; the shell wrapper loops):
#
#   Rscript 06_260526-drimseq_pseudotime_pseudobulk.R \
#       --input  /path/to/merged_*.h5ad \
#       --lineage inhibitory \
#       --pseudotime_col dpt_pseudotime_inhibitory \
#       --output_dir /path/to/results/inhibitory/pseudobulk_pseudotime \
#       --keep_induced_false_only \
#       --n_bins 3 --min_cells_per_unit 30 --ncores 16
#
# Test mode:
#   Rscript 06_260526-drimseq_pseudotime_pseudobulk.R --test
# =============================================================================

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT") # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
  library(rlang)
  library(Matrix)
  library(splines)
  library(jsonlite) # run_metadata.json sidecar
})

reticulate::use_condaenv("apa-r", required = TRUE)
box::use(
  DRIMSeq[
    dmDSdata, dmFilter, dmPrecision, dmFit, dmTest, results,
    samples, counts, common_precision
  ],
  dplyr[
    arrange, mutate, select, filter, `%>%`, bind_rows, left_join,
    group_by, summarise, n_distinct, n
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
    pseudotime_col = NULL, # default derived from --lineage if NULL
    # Batch covariate for the fixed-effect adjustment.
    # Default: batch_number (= sample_id = "library" in Li et al. 2023 Nature,
    # which uses `library` as the per-sample covariate in their DE model:
    # Y ~ n_UMI + library + cell_type + condition).
    # The coarser super_batch (3 levels; sequencing run) is not used
    # for the controls-only pseudotime test: in controls super_batch has
    # only 2 levels and the `213` level contains a single sample (213664),
    # so super_batch fails to control for between-sample biological
    # variation among the 7 samples in super_batch 206. The cell-level LRT
    # and the permutation test also use batch_number as the
    # covariate.
    covariate_col = "batch_number",
    output_dir = NULL,
    n_bins = 3L, # 3 bins (vs 8 in the treatment script)
    min_cells_per_unit = 30L, # raise from 10 (3-bin design has 10x more cells/unit)
    min_samps_gene_expr = 5L,
    min_gene_expr = 50L,
    min_samps_feature_expr = 5L,
    min_feature_expr = 20L,
    min_units_per_bin = 4L, # require >= 4 sample-bin units per bin with >= min_bin_expr counts
    min_bin_expr = 10L,
    keep_induced_false_only = FALSE, # set TRUE for Knoblich controls-only design
    seed = 42L,
    ncores = max(1L, parallel::detectCores() - 1L),
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
    if (a == "--keep_induced_false_only") {
      defaults$keep_induced_false_only <- TRUE
      i <- i + 1L
      next
    }
    if (startsWith(a, "--") && i < length(argv)) {
      key <- sub("^--", "", a)
      val <- argv[i + 1L]
      if (!(key %in% names(defaults))) stop("Unknown argument: ", a)
      cur <- defaults[[key]]
      defaults[[key]] <- if (!is.null(cur) && is.integer(cur)) as.integer(val) else val
      i <- i + 2L
    } else {
      stop("Bad argument: ", a)
    }
  }
  if (!defaults$test && is.null(defaults$output_dir)) {
    stop("--output_dir is required when not in --test mode")
  }
  # Default pseudotime_col: derive from --lineage if not explicitly set.
  # Special case: --lineage all (mouse single-trajectory convention)
  # falls back to the bare "dpt_pseudotime" column, not "dpt_pseudotime_all"
  # (which does not exist in either dataset).
  if (is.null(defaults$pseudotime_col)) {
    if (defaults$lineage %in% c("all", "")) {
      defaults$pseudotime_col <- "dpt_pseudotime"
    } else {
      defaults$pseudotime_col <- paste0("dpt_pseudotime_", defaults$lineage)
    }
  }
  if (defaults$n_bins < 2L) {
    stop("--n_bins must be >= 2 (need at least 2 bins for the LRT to have df > 0)")
  }
  defaults
}

# =============================================================================
# Section 1: Helper functions
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

# Note: the treatment-arm pseudobulk (02_drimseq_treatment_main_pseudobulk.R) has
# a `derive_super_batch <- function(sample_id_vec) substr(sample_id_vec, 1, 3)`
# helper because there `induced` is collinear with sample, forcing the coarser
# `super_batch` (sequencing run) as the only viable covariate. The pseudotime
# script (this file) does not need that helper: the covariate is
# `batch_number` (= sample_id = Li et al. 2023 "library"), and `super_batch`
# would be a degenerate covariate in the controls-only subset (see the
# comment on covariate_col in parse_args()). A treatment-arm variant of this
# script would need the helper again.

to_dgCMatrix_double <- function(X) {
  if (inherits(X, "RsparseMatrix")) {
    tmp <- as(X, "TsparseMatrix")
    return(Matrix::sparseMatrix(
      i = tmp@i + 1L, j = tmp@j + 1L,
      x = as.numeric(tmp@x),
      dims = dim(tmp), dimnames = dimnames(tmp)
    ))
  }
  if (!inherits(X, "dgCMatrix")) X <- as(X, "CsparseMatrix")
  if (typeof(X@x) != "double") X@x <- as.numeric(X@x)
  X
}

# For the pseudotime design (no binary contrast), require >= min_units_per_bin
# sample-bin units per pt_bin to have gene-level counts >= min_bin_expr.
# Prevents single-bin-dominated genes (e.g. a gene that is only expressed in
# bin 0 has no leverage for a pt_bin effect test).
min_units_per_bin_prefilter <- function(counts_df, samples_df,
                                        min_units_per_bin, min_bin_expr) {
  sample_cols <- samples_df$sample_id
  stopifnot(all(sample_cols %in% colnames(counts_df)))
  bin_ids <- sort(unique(samples_df$bin_id))
  if (length(bin_ids) < 2L) {
    stop("min_units_per_bin_prefilter: need >= 2 pt_bin levels.")
  }
  keep_genes <- rep(TRUE, length(unique(counts_df$gene_id)))
  n_units_by_bin <- list()
  for (b in bin_ids) {
    cols_b <- sample_cols[samples_df$bin_id == b]
    if (length(cols_b) == 0L) {
      keep_genes <- rep(FALSE, length(unique(counts_df$gene_id)))
      names(keep_genes) <- unique(counts_df$gene_id)
      n_units_by_bin[[as.character(b)]] <- 0
      next
    }
    ct_b <- as.matrix(counts_df[, cols_b, drop = FALSE])
    storage.mode(ct_b) <- "double"
    gid <- counts_df$gene_id
    gene_b <- rowsum(ct_b, group = gid, reorder = FALSE)
    n_units_b <- rowSums(gene_b >= min_bin_expr)
    n_units_by_bin[[as.character(b)]] <- n_units_b
  }
  # Per-gene min over bins of n_units >= min_units_per_bin
  n_units_mat <- do.call(cbind, n_units_by_bin)
  min_n_units <- apply(n_units_mat, 1L, min)
  keep_genes <- min_n_units >= min_units_per_bin
  kept_gene_ids <- rownames(n_units_mat)[keep_genes]
  keep_rows <- counts_df$gene_id %in% kept_gene_ids
  list(
    counts_df_kept = counts_df[keep_rows, , drop = FALSE],
    n_genes_pre = nrow(n_units_mat), n_genes_post = sum(keep_genes),
    n_features_pre = nrow(counts_df), n_features_post = sum(keep_rows)
  )
}

quantile_bin_pseudotime <- function(pt, n_bins) {
  qs <- stats::quantile(pt,
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE, names = FALSE
  )
  qs <- unique(qs)
  if (length(qs) < 2L) stop("Pseudotime has fewer than 2 unique values; cannot bin.")
  bin <- cut(pt, breaks = qs, include.lowest = TRUE, labels = FALSE)
  bin_mid <- vapply(seq_len(max(bin, na.rm = TRUE)), function(b) {
    mean(pt[bin == b], na.rm = TRUE)
  }, numeric(1))
  list(bin = bin, bin_mid = bin_mid)
}

# Pseudotime design: full = ~ <covariate> + factor(pt_bin); null = ~ <covariate>.
# LRT df per (K-1) compositional dimension = n_bins - 1. With n_bins = 3 this
# gives a 2-df test that can detect any 3-bin pattern (monotonic, U, inverted-U).
# `spline_df` is unused (kept for API parity with the treatment script).
# `covariate_col` is required (no default). Callers must pass it explicitly;
# this prevents silent fallback to a stale default when the script is invoked
# from a shell driver that intends a different covariate. The Knoblich and
# mouse drivers both pass --covariate_col batch_number; the test suite (Tests
# 3 and 6) also exercises the production path with batch_number.
build_design_matrices <- function(samples_table, spline_df = NULL,
                                  covariate_col) {
  if (missing(covariate_col) || is.null(covariate_col)) {
    stop("build_design_matrices: covariate_col is required.")
  }
  if (!(covariate_col %in% colnames(samples_table))) {
    stop(
      "build_design_matrices: covariate_col '", covariate_col,
      "' not found in samples_table. Available: ",
      paste(colnames(samples_table), collapse = ", ")
    )
  }
  s <- samples_table
  s[[covariate_col]] <- factor(as.character(s[[covariate_col]]))
  s$pt_bin <- factor(as.character(s$bin_id))
  fml_full <- as.formula(paste0("~ ", covariate_col, " + pt_bin"))
  fml_null <- as.formula(paste0("~ ", covariate_col))
  list(
    full = model.matrix(fml_full, data = s),
    null = model.matrix(fml_null, data = s),
    covariate_col = covariate_col
  )
}

# Delta_psi pseudotime: signed Wasserstein-1 between bin n_bins-1 and bin 0
# mean proportions. Sign convention: positive = lengthening (distal shift)
# from early bin to late bin. Same rank-normalized W as the treatment-arm
# delta_psi; just contrasts the extreme bins (bin0 vs bin_last) instead of
# induced True vs False. Per-bin psi values are also reported as
# psi_bin0..psi_bin{n_bins-1} for trajectory plotting.
compute_delta_psi_pseudotime <- function(counts_d, samples_df) {
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

  bin_levels <- sort(unique(samples_df$bin_id))
  n_bins_eff <- length(bin_levels)
  bin_first <- bin_levels[1L]
  bin_last <- bin_levels[n_bins_eff]
  bin_str <- function(b) paste0("psi_bin", b)
  n_cells_str <- function(b) paste0("n_cells_bin", b)
  n_units_str <- function(b) paste0("n_units_bin", b)

  psi_per_bin <- vector("list", n_bins_eff)
  n_cells_per_bin <- vector("list", n_bins_eff)
  n_units_per_bin <- vector("list", n_bins_eff)
  feat_pi_first <- numeric(nrow(feature_meta))
  feat_pi_last <- numeric(nrow(feature_meta))

  n_cells <- as.numeric(samples_df$n_cells)
  bin_id_v <- samples_df$bin_id

  for (bi in seq_along(bin_levels)) {
    psi_per_bin[[bi]] <- rep(NA_real_, length(gene_unique))
    n_cells_per_bin[[bi]] <- rep(NA_real_, length(gene_unique))
    n_units_per_bin[[bi]] <- rep(0L, length(gene_unique))
  }

  for (gi in seq_along(gene_unique)) {
    g <- gene_unique[gi]
    rows <- which(gene_id == g)
    K <- length(rows)
    if (K < 2L) {
      for (bi in seq_along(bin_levels)) {
        psi_per_bin[[bi]][gi] <- NA_real_
        n_cells_per_bin[[bi]][gi] <- NA_real_
        n_units_per_bin[[bi]][gi] <- 0L
      }
      feat_pi_first[rows] <- NA_real_
      feat_pi_last[rows] <- NA_real_
      next
    }
    rank_norm <- (seq_len(K) - 1) / (K - 1)
    Y_g <- counts_mat[rows, , drop = FALSE]
    unit_totals <- colSums(Y_g)
    keep_u <- unit_totals > 0
    if (sum(keep_u) == 0L) {
      for (bi in seq_along(bin_levels)) {
        psi_per_bin[[bi]][gi] <- NA_real_
        n_cells_per_bin[[bi]][gi] <- NA_real_
        n_units_per_bin[[bi]][gi] <- 0L
      }
      feat_pi_first[rows] <- NA_real_
      feat_pi_last[rows] <- NA_real_
      next
    }
    Y_keep <- Y_g[, keep_u, drop = FALSE]
    tot_keep <- unit_totals[keep_u]
    prop_uk <- sweep(Y_keep, 2L, tot_keep, "/")
    psi_u <- as.numeric(crossprod(rank_norm, prop_uk))
    bin_id_keep <- bin_id_v[keep_u]
    w_keep <- n_cells[keep_u]
    for (bi in seq_along(bin_levels)) {
      b <- bin_levels[bi]
      in_b <- bin_id_keep == b
      w_b <- w_keep * in_b
      sum_w_b <- sum(w_b)
      if (sum_w_b > 0) {
        psi_per_bin[[bi]][gi] <- sum(psi_u * w_b) / sum_w_b
        if (b == bin_first) {
          feat_pi_first[rows] <- as.numeric(prop_uk %*% (w_b / sum_w_b))
        }
        if (b == bin_last) {
          feat_pi_last[rows] <- as.numeric(prop_uk %*% (w_b / sum_w_b))
        }
      } else {
        psi_per_bin[[bi]][gi] <- NA_real_
        if (b == bin_first) feat_pi_first[rows] <- NA_real_
        if (b == bin_last) feat_pi_last[rows] <- NA_real_
      }
      n_cells_per_bin[[bi]][gi] <- sum_w_b
      n_units_per_bin[[bi]][gi] <- sum(in_b)
    }
  }
  # Assemble gene-level tibble
  gene_tbl <- tibble::tibble(gene_id = gene_unique)
  for (bi in seq_along(bin_levels)) {
    b <- bin_levels[bi]
    gene_tbl[[bin_str(b)]] <- psi_per_bin[[bi]]
    gene_tbl[[n_cells_str(b)]] <- n_cells_per_bin[[bi]]
    gene_tbl[[n_units_str(b)]] <- n_units_per_bin[[bi]]
  }
  psi_first <- psi_per_bin[[1L]]
  psi_last <- psi_per_bin[[n_bins_eff]]
  gene_tbl$delta_psi_pseudotime <- psi_last - psi_first
  # Feature-level pi (only first vs last bin to keep the table small)
  feature_tbl <- tibble::tibble(
    gene_id = feature_meta$gene_id,
    feature_id = feature_meta$feature_id,
    pi_bin_first = feat_pi_first,
    pi_bin_last = feat_pi_last,
    delta_pi_feature = feat_pi_last - feat_pi_first
  )
  list(gene = gene_tbl, feature = feature_tbl)
}

# =============================================================================
# Section 2: DRIMSeq pipeline wrapper: single LRT, gene-wise precision with
# NaN repair.
# =============================================================================
# Precision strategy: leave gene-wise precision as DRIMSeq estimates it, but
# replace NaN/Inf/<=0 entries (typically K>=7 in this dataset) with the
# common precision. This light repair suffices for designs where gene-wise
# estimates are well-behaved for most genes (verified 2026-05-13 for the
# treatment main-effect design on full inhibitory: 5,291/5,441
# valid p-values, 76 hits at adj_p<0.01).
run_drimseq_pipeline <- function(d, design_full, design_null,
                                 bpp, add_uniform = TRUE,
                                 seed = 42L, verbose = TRUE) {
  if (verbose) message("  dmPrecision (genewise, prec_subset=1, seed=", seed, ")...")
  used_genewise <- TRUE
  # Pin prec_subset = 1 for deterministic precision estimation
  # (dmPrecision is non-deterministic by default; common_precision
  # varied 95-984 in repeated runs of the same data). Seed before each call
  # because add_uniform=TRUE adds tiny jitter to zero-count features.
  set.seed(seed)
  d <- tryCatch(
    DRIMSeq::dmPrecision(d,
      design = design_full,
      prec_subset = 1,
      add_uniform = add_uniform, BPPARAM = bpp
    ),
    error = function(e) {
      if (verbose) {
        message("  genewise precision failed: ", conditionMessage(e))
        message("  retrying with common precision only...")
      }
      used_genewise <<- FALSE
      set.seed(seed)
      DRIMSeq::dmPrecision(d,
        design = design_full,
        prec_subset = 1,
        genewise_precision = FALSE,
        add_uniform = add_uniform, BPPARAM = bpp
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
        if (verbose) message("  precision: dmPrecision already ran with genewise=FALSE; skipping NaN repair")
        "skipped_common"
      } else {
        gw_vec <- d@genewise_precision
        cp <- as.numeric(d@common_precision)[1]
        stopifnot(
          is.numeric(gw_vec), length(gw_vec) > 0L,
          is.finite(cp), cp > 0
        )
        bad <- !is.finite(gw_vec) | gw_vec <= 0
        n_repaired <<- as.integer(sum(bad))
        if (n_repaired > 0L) {
          gw_new <- gw_vec
          gw_new[bad] <- cp
          d@genewise_precision <- gw_new
          if (verbose) {
            message(
              "  precision NaN repair: ", n_repaired, "/",
              length(gw_vec),
              " gene-wise precisions had NaN/Inf/<=0; ",
              "replaced with common_precision (",
              round(cp, 4), ")"
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
              length(gw_vec), " have finite positive precision)"
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

  if (verbose) message("  dmFit...")
  d <- DRIMSeq::dmFit(d, design = design_full, BPPARAM = bpp)

  if (verbose) message("  dmTest (treatment main)...")
  d_test <- DRIMSeq::dmTest(d, design = design_null, BPPARAM = bpp)

  list(
    d_fit = d,
    gene = DRIMSeq::results(d_test, level = "gene") |>
      tibble::as_tibble() |> arrange(pvalue),
    feature = DRIMSeq::results(d_test, level = "feature") |>
      tibble::as_tibble() |> arrange(pvalue),
    common_precision = tryCatch(
      as.numeric(DRIMSeq::common_precision(d)),
      error = function(e) NA_real_
    ),
    n_genes_precision_repaired = n_repaired
  )
}

# =============================================================================
# Section 3: no cluster bootstrap in the pseudotime script.
# =============================================================================
# The treatment-arm script keeps a cluster_bootstrap function as a directional
# sensitivity analysis for the induced contrast (where the effective n_samples
# is 22, borderline for cluster bootstrap reliability).
#
# The pseudotime test does not need it: there is no binary contrast, the n_eff
# is the number of pseudobulk units (24 per lineage), and the LRT is testing
# a multi-level factor whose stability under sample bootstrap would just
# repeat what the permutation test already establishes. As a
# resampling-based sensitivity check, the within-sample permutation test
# (04_260526-permutation_optimized_2.R) is the more powerful
# alternative.
# =============================================================================
# (The cluster-bootstrap implementation is in
# ../perturbation/02_drimseq_treatment_main_pseudobulk.R.)

# =============================================================================
# Section 4: Test mode (--test)
# =============================================================================

run_tests <- function() {
  cat("=== drimseq_pseudotime_pseudobulk: verification tests ===\n\n")
  all_pass <- TRUE

  cat("Test 1: pseudobulk_by_group correctness... ")
  set.seed(42)
  n_cells <- 100L
  n_feat <- 7L
  dense_in <- matrix(as.numeric(rpois(n_cells * n_feat, 3)), nrow = n_cells)
  cells_x_feat <- as(Matrix::Matrix(dense_in, sparse = TRUE), "CsparseMatrix")
  grp <- sample(letters[1:4], n_cells, replace = TRUE)
  pb <- pseudobulk_by_group(cells_x_feat, grp)
  expected <- vapply(sort(unique(grp)), function(g) {
    colSums(as.matrix(cells_x_feat[grp == g, , drop = FALSE]))
  }, numeric(n_feat))
  rownames(expected) <- NULL
  if (max(abs(as.matrix(pb) - expected)) < 1e-10) {
    cat("PASS\n")
  } else {
    cat("FAIL\n")
    all_pass <- FALSE
  }

  # (No test of derive_super_batch, which is not part of this script.
  # super_batch is a treatment-arm concept; the pseudotime script uses
  # batch_number unconditionally.)

  cat("Test 2: quantile_bin_pseudotime (3 bins)... ")
  bres <- quantile_bin_pseudotime(runif(900), 3L)
  tab <- table(bres$bin)
  if (length(tab) == 3L && max(tab) - min(tab) <= 2L) {
    cat("PASS (sizes:", paste(as.integer(tab), collapse = ","), ")\n")
  } else {
    cat("FAIL (sizes:", paste(as.integer(tab), collapse = ","), ")\n")
    all_pass <- FALSE
  }

  cat("Test 3: design-matrix df arithmetic (pseudotime, 3 bins, batch_number covariate)... ")
  n_units <- 24L # 8 samples x 3 bins, as in the Knoblich controls
  samples_table <- data.frame(
    sample_id = sprintf("u%03d", seq_len(n_units)),
    bin_id = rep(1L:3L, times = 8L),
    # batch_number = sample_id 1-to-1 (the production case): 8 batches
    batch_number = rep(sprintf("B%02d", seq_len(8L)), each = 3L),
    stringsAsFactors = FALSE
  )
  designs <- build_design_matrices(samples_table, covariate_col = "batch_number")
  # full: 1 (int) + 7 (batch-1) + 2 (pt_bin-1) = 10
  # null: 1 (int) + 7 (batch-1)               = 8
  # df per (K-1) compositional dimension = 2
  df_main <- ncol(designs$full) - ncol(designs$null)
  if (ncol(designs$full) == 10L && ncol(designs$null) == 8L && df_main == 2L) {
    cat("PASS (full=", ncol(designs$full), " null=", ncol(designs$null),
      " df=", df_main, ")\n",
      sep = ""
    )
  } else {
    cat("FAIL (full=", ncol(designs$full), " null=", ncol(designs$null),
      " df=", df_main, ")\n",
      sep = ""
    )
    all_pass <- FALSE
  }

  cat("Test 4: compute_delta_psi_pseudotime sign convention... ")
  # Gene A: bin 1 has all UTR_first, bin 3 has all UTR_third (= distal)
  # rank_norm for 3 features = (0, 0.5, 1)
  # psi_bin1 = 1*0 + 0*0.5 + 0*1 = 0
  # psi_bin3 = 0*0 + 0*0.5 + 1*1 = 1
  # delta_psi_pseudotime = psi_bin3 - psi_bin1 = +1 (full lengthening)
  counts_d <- data.frame(
    gene_id = c("A", "A", "A", "B", "B"),
    feature_id = c("A1", "A2", "A3", "B1", "B2"),
    u1 = c(80, 10, 10, 50, 50),
    u2 = c(40, 30, 30, 50, 50),
    u3 = c(10, 10, 80, 50, 50),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  samples_df <- data.frame(
    sample_id = c("u1", "u2", "u3"),
    bin_id = c(1L, 2L, 3L),
    n_cells = c(100, 100, 100),
    stringsAsFactors = FALSE
  )
  dp <- compute_delta_psi_pseudotime(counts_d, samples_df)
  a_row <- dp$gene[dp$gene$gene_id == "A", ]
  expected_dpsi <- 0.7 # psi_bin3=0.7, psi_bin1=0.1 (since 0.1*0 + 0.1*0.5 + 0.8*1=0.85; etc)
  # Quick analytical check: bin1 fracs = (80,10,10)/100; psi = 0*0.8 + 0.5*0.1 + 1*0.1 = 0.15
  # bin3 fracs = (10,10,80)/100; psi = 0*0.1 + 0.5*0.1 + 1*0.8 = 0.85
  # delta_psi = 0.85 - 0.15 = 0.70
  if (!is.na(a_row$delta_psi_pseudotime) &&
    abs(a_row$delta_psi_pseudotime - 0.70) < 1e-6) {
    cat("PASS (delta_psi =", round(a_row$delta_psi_pseudotime, 4), ")\n")
  } else {
    cat("FAIL (delta_psi =", a_row$delta_psi_pseudotime, ")\n")
    all_pass <- FALSE
  }

  cat("Test 5: min_units_per_bin_prefilter drops single-bin-only genes... ")
  counts_d6 <- data.frame(
    gene_id = c(rep("A", 2L), rep("B", 2L)),
    feature_id = c("A1", "A2", "B1", "B2"),
    u1 = c(100, 50, 100, 50),
    u2 = c(100, 50, 0, 0),
    u3 = c(100, 50, 0, 0),
    u4 = c(100, 50, 0, 0),
    u5 = c(100, 50, 0, 0),
    u6 = c(100, 50, 0, 0),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  samples_df6 <- data.frame(
    sample_id = c("u1", "u2", "u3", "u4", "u5", "u6"),
    bin_id = c(1L, 2L, 3L, 1L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  pre <- min_units_per_bin_prefilter(counts_d6, samples_df6,
    min_units_per_bin = 1L, min_bin_expr = 10L
  )
  # Gene A has expression in all bins; Gene B has expression only in bin 1
  if (pre$n_genes_post == 1L && pre$counts_df_kept$gene_id[1] == "A") {
    cat("PASS (kept", pre$n_genes_post, "of", pre$n_genes_pre, "genes)\n")
  } else {
    cat(
      "FAIL (kept", pre$n_genes_post, "of", pre$n_genes_pre, "genes; kept IDs: ",
      paste(unique(pre$counts_df_kept$gene_id), collapse = ","), ")\n"
    )
    all_pass <- FALSE
  }

  cat("Test 6: end-to-end pipeline runs under H0 (no pt_bin effect, batch_number covariate)... ")
  set.seed(42L)
  n_genes <- 100L
  n_bins_t <- 3L
  pi_true <- c(0.3, 0.3, 0.4)
  K_t <- length(pi_true)
  s_true <- 50
  # 8 samples; batch_number 1-to-1 with sample (the production setup; matches
  # Knoblich 8 controls and is conservative for mouse's 6).
  sample_meta <- data.frame(
    sample_id_root = sprintf("S%02d", seq_len(8L)),
    batch_number = sprintf("B%02d", seq_len(8L)),
    stringsAsFactors = FALSE
  )
  units <- expand.grid(
    sample_id_root = sample_meta$sample_id_root,
    bin_id = seq_len(n_bins_t),
    stringsAsFactors = FALSE
  )
  units$sample_id <- paste(units$sample_id_root, units$bin_id, sep = "__")
  units <- merge(units, sample_meta, by = "sample_id_root")
  units$n_cells <- 100L
  n_units_t <- nrow(units)
  counts_list <- vector("list", n_genes)
  alpha <- pi_true * s_true
  for (gi in seq_len(n_genes)) {
    Y <- t(vapply(seq_len(n_units_t), function(i) {
      p_dir <- rgamma(K_t, alpha, 1)
      p_dir <- p_dir / sum(p_dir)
      as.integer(rmultinom(1, units$n_cells[i], p_dir))
    }, integer(K_t)))
    df_g <- data.frame(
      gene_id = paste0("G", gi),
      feature_id = paste0("G", gi, "_", seq_len(K_t)),
      t(Y), check.names = FALSE, stringsAsFactors = FALSE
    )
    colnames(df_g)[3:(2 + n_units_t)] <- units$sample_id
    counts_list[[gi]] <- df_g
  }
  counts_df <- do.call(rbind, counts_list)
  samples_drim <- data.frame(
    sample_id = units$sample_id, sample_id_root = units$sample_id_root,
    bin_id = units$bin_id, batch_number = units$batch_number,
    n_cells = units$n_cells, stringsAsFactors = FALSE
  )
  bpp <- BiocParallel::SerialParam()
  d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples_drim)
  d <- DRIMSeq::dmFilter(d,
    min_samps_gene_expr = 5L, min_gene_expr = 50L,
    min_samps_feature_expr = 5L, min_feature_expr = 20L
  )
  designs <- build_design_matrices(DRIMSeq::samples(d), covariate_col = "batch_number")
  res <- run_drimseq_pipeline(d, designs$full, designs$null,
    bpp = bpp, seed = 42L, verbose = FALSE
  )
  pvals <- res$gene$pvalue[!is.na(res$gene$pvalue)]
  mean_p <- if (length(pvals) > 0) mean(pvals) else NA_real_
  if (length(pvals) / n_genes >= 0.5 && !is.na(mean_p) &&
    mean_p >= 0.25 && mean_p <= 0.75) {
    cat("PASS (n_valid=", length(pvals), " mean_p=", round(mean_p, 3), ")\n", sep = "")
  } else {
    cat("FAIL (n_valid=", length(pvals), " mean_p=", round(mean_p, 3), ")\n", sep = "")
    all_pass <- FALSE
  }

  cat("Test 7: reproducibility under fixed seed... ")
  res_a <- run_drimseq_pipeline(d, designs$full, designs$null,
    bpp = bpp, seed = 42L, verbose = FALSE
  )
  res_b <- run_drimseq_pipeline(d, designs$full, designs$null,
    bpp = bpp, seed = 42L, verbose = FALSE
  )
  if (identical(res_a$gene$pvalue, res_b$gene$pvalue) &&
    identical(res_a$common_precision, res_b$common_precision)) {
    cat(
      "PASS (",
      length(res_a$gene$pvalue),
      "p-values identical, common_precision =",
      round(res_a$common_precision, 4), ")\n"
    )
  } else {
    cat("FAIL\n")
    all_pass <- FALSE
  }

  cat("\n=== ", if (all_pass) "ALL TESTS PASS" else "TESTS FAILED", " ===\n", sep = "")
  if (!all_pass) quit(status = 1L)
}

# =============================================================================
# Main
# =============================================================================

opts <- parse_args(commandArgs(trailingOnly = TRUE))
if (opts$test) {
  run_tests()
  quit(status = 0L)
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
run_start_time <- format(Sys.time())
run_host <- Sys.info()[["nodename"]]
bpp <- BiocParallel::MulticoreParam(workers = opts$ncores)
slug <- gsub("[^A-Za-z0-9]+", "_", tolower(opts$lineage))

cat("=== DRIMSeq pseudotime LRT (pseudobulk, 3-bin factor; 260526) ===\n")
cat("  input:                ", opts$input, "\n")
cat("  output_dir:           ", opts$output_dir, "\n")
cat("  lineage:              ", opts$lineage, "\n")
cat("  pseudotime_col:       ", opts$pseudotime_col, "\n")
cat("  n_bins:               ", opts$n_bins, "\n")
cat("  min_cells_per_unit:   ", opts$min_cells_per_unit, "\n")
cat("  min_units_per_bin:    ", opts$min_units_per_bin, "\n")
cat("  min_bin_expr:         ", opts$min_bin_expr, "\n")
cat("  keep_induced_false_only:", opts$keep_induced_false_only, "\n")
cat("  ncores:               ", opts$ncores, "\n")
cat("  seed:                 ", opts$seed, "\n")
cat("  start:                ", run_start_time, "\n")

# ---- Load + lineage subset + controls filter ----
cat("\n[Load] Reading h5ad...\n")
adata <- anndata::read_h5ad(opts$input)
required_obs <- c("sample_id", opts$pseudotime_col)
if (opts$keep_induced_false_only) {
  required_obs <- c(required_obs, "induced")
}
missing_obs <- setdiff(required_obs, colnames(adata$obs))
if (length(missing_obs) > 0L) stop("Missing obs columns: ", paste(missing_obs, collapse = ", "))
if (!"gene_name" %in% colnames(adata$var)) stop("Missing var column gene_name")

obs_all <- adata$obs
var_all <- adata$var
pt_all <- as.numeric(obs_all[[opts$pseudotime_col]])
sample_id_all <- as.character(obs_all$sample_id)
# Covariate source: must be a column already present in obs. (There is no
# substr-of-sample_id "super_batch" derivation as in the treatment-arm
# script; the covariate here is batch_number.)
if (!(opts$covariate_col %in% colnames(obs_all))) {
  stop(
    "Covariate column '", opts$covariate_col,
    "' not found in obs; available columns: ",
    paste(colnames(obs_all), collapse = ", ")
  )
}
covariate_all <- as.character(obs_all[[opts$covariate_col]])
cat("  covariate '", opts$covariate_col, "' read from obs (",
  length(unique(covariate_all)), " unique levels)\n",
  sep = ""
)
induced_all <- if ("induced" %in% colnames(obs_all)) {
  as.character(obs_all$induced)
} else {
  rep(NA_character_, length(pt_all))
}
keep <- !is.na(pt_all)
cat("  cells in lineage (non-NA pt):", sum(keep), "/", length(keep), "\n")

if (opts$keep_induced_false_only) {
  is_false <- induced_all %in% c("FALSE", "False", "false", "F", "0")
  keep <- keep & is_false
  cat("  after keep_induced_false_only:", sum(keep), "\n")
}

X_all <- to_dgCMatrix_double(adata$X)
X_lin <- X_all[keep, , drop = FALSE]
pt_lin <- pt_all[keep]
sample_id_lin <- sample_id_all[keep]
covariate_lin <- covariate_all[keep]
rm(X_all, adata)
gc(verbose = FALSE)

sample_table_lin <- data.frame(
  sample_id = sample_id_lin,
  stringsAsFactors = FALSE
)
sample_table_lin[[opts$covariate_col]] <- covariate_lin
cat("\nSample composition (kept cells):\n")
print(table(unique(sample_table_lin)[[opts$covariate_col]],
  dnn = opts$covariate_col
))
cat("  samples (n =", length(unique(sample_id_lin)), "):\n")
print(table(sample_id_lin))

# ---- Pseudobulk ----
cat("\n[Pseudobulk] Quantile-binning pseudotime into", opts$n_bins, "bins (global)...\n")
set.seed(opts$seed) # belt-and-braces, no randomness in quantile_bin_pseudotime
bres <- quantile_bin_pseudotime(pt_lin, opts$n_bins)
n_bins_obs <- max(bres$bin, na.rm = TRUE)
if (n_bins_obs < opts$n_bins) {
  cat(
    "  WARN: only", n_bins_obs, "distinct bin levels produced (asked for",
    opts$n_bins, "); proceeding with", n_bins_obs, "bins.\n"
  )
}
cat("  bin sizes:", paste(as.integer(table(bres$bin)), collapse = ", "), "\n")
cat("  bin midpoints:", paste(round(bres$bin_mid, 3), collapse = ", "), "\n")

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
pb_meta <- left_join(pb_meta, sample_lookup,
  by = c("sample_id_root" = "sample_id")
)
pb_meta$n_cells <- as.integer(table(groups)[pb_meta$sample_id])
keep_units <- pb_meta$n_cells >= opts$min_cells_per_unit
cat("  units after min_cells_per_unit (>=", opts$min_cells_per_unit, "): ",
  sum(keep_units), "/", nrow(pb_meta), "\n",
  sep = ""
)
pb_counts <- pb_counts[, keep_units, drop = FALSE]
pb_meta <- pb_meta[keep_units, , drop = FALSE]

manifest_path <- file.path(opts$output_dir, "pseudobulk_units.csv")
write.csv(as.data.frame(pb_meta), manifest_path, row.names = FALSE)
# Diagnostic and DRIMSeq-input column selection must use the actual covariate
# column requested via --covariate_col (default batch_number), not a hardcoded
# "super_batch" as in the treatment-arm pseudobulk
# script where super_batch is the only viable covariate (the `induced`
# contrast is collinear with sample). In the pseudotime script that would be wrong
# whenever --covariate_col != "super_batch" (i.e. essentially always for
# Knoblich controls and always for mouse), because pb_meta only carries the
# requested covariate column. The stopifnot below checks
# that the covariate is in pb_meta; if it fails, the upstream
# left_join from sample_lookup is broken.
stopifnot(opts$covariate_col %in% colnames(pb_meta))
cat("  manifest:", manifest_path, "\n  unit composition (",
  opts$covariate_col, " x bin):\n",
  sep = ""
)
print(table(pb_meta[[opts$covariate_col]], pb_meta$bin_id,
  dnn = c(opts$covariate_col, "bin_id")
))

# ---- DRIMSeq pipeline ----
feature_ids <- rownames(var_all)
stopifnot(!anyDuplicated(feature_ids))
counts_df <- data.frame(
  gene_id = as.character(var_all$gene_name),
  feature_id = feature_ids, pb_counts,
  check.names = FALSE, stringsAsFactors = FALSE
)
pb_samples_drim <- as.data.frame(pb_meta)[, c(
  "sample_id", "sample_id_root", "bin_id", "pt_mid",
  opts$covariate_col, "n_cells"
)]
stopifnot(all(colnames(pb_counts) == pb_samples_drim$sample_id))

dm_filter_args <- list(
  min_samps_gene_expr    = opts$min_samps_gene_expr,
  min_gene_expr          = opts$min_gene_expr,
  min_samps_feature_expr = opts$min_samps_feature_expr,
  min_feature_expr       = opts$min_feature_expr
)

cat("\n[Pre-filter] min_units_per_bin pre-filter (>=",
  opts$min_units_per_bin, "units per bin with >=", opts$min_bin_expr,
  "counts)...\n",
  sep = ""
)
prefilt <- min_units_per_bin_prefilter(
  counts_df, pb_samples_drim,
  opts$min_units_per_bin, opts$min_bin_expr
)
cat("  pre :", prefilt$n_genes_pre, "genes,", prefilt$n_features_pre, "features\n")
cat("  post:", prefilt$n_genes_post, "genes,", prefilt$n_features_post, "features\n")
if (prefilt$n_features_post == 0L) stop("Pre-filter dropped all genes")
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
if (n_feat_post == 0L) stop("dmFilter dropped all features")

cat("\n[DRIMSeq] Building pseudotime design (~ ", opts$covariate_col,
  " + factor(pt_bin))...\n",
  sep = ""
)
designs <- build_design_matrices(DRIMSeq::samples(d),
  covariate_col = opts$covariate_col
)
cat("  full design:", ncol(designs$full), "cols (rank=", qr(designs$full)$rank, ")\n", sep = "")
cat("  null design:", ncol(designs$null), "cols\n")
cat("  LRT df per (K-1):", ncol(designs$full) - ncol(designs$null), "\n")

cat("\n[DRIMSeq] Fit + 1 test (pseudotime main, no induced contrast) on observed data...\n")
obs_out <- run_drimseq_pipeline(d, designs$full, designs$null,
  bpp = bpp, seed = opts$seed, verbose = TRUE
)
cat("  common precision =", round(obs_out$common_precision, 4), "\n")
cat("  precision NaN-repaired:", obs_out$n_genes_precision_repaired, "/",
  n_genes_post, "\n",
  sep = ""
)

# ---- Delta_psi (pseudotime: bin_last - bin_first) ----
cat("\n[Delta_psi] Computing per-bin psi + bin0-vs-bin{last} delta_psi_pseudotime...\n")
dpsi <- compute_delta_psi_pseudotime(counts_post, DRIMSeq::samples(d))
dpsi_gene_path <- file.path(opts$output_dir, sprintf("delta_psi_pseudotime_gene_%s.csv", slug))
dpsi_feat_path <- file.path(opts$output_dir, sprintf("delta_psi_pseudotime_feature_%s.csv", slug))
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
  mutate(
    lineage = opts$lineage,
    test_type = "pseudotime_main",
    master_seed = as.integer(opts$seed),
    n_bins = as.integer(opts$n_bins),
    min_cells_per_unit = as.integer(opts$min_cells_per_unit),
    covariate_col = opts$covariate_col
  )
feat_df <- obs_out$feature %>%
  left_join(dpsi$feature, by = c("gene_id", "feature_id")) %>%
  mutate(lineage = opts$lineage, test_type = "pseudotime_main")

gene_csv <- file.path(opts$output_dir, sprintf("drimseq_gene_%s_pseudotime.csv", slug))
feat_csv <- file.path(opts$output_dir, sprintf("drimseq_feature_diagnostic_%s_pseudotime.csv", slug))
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

# ---- Summary CSV ----
summary_row <- tibble(
  lineage = opts$lineage, n_cells_lineage = sum(keep),
  n_units_pseudobulk = nrow(pb_samples_drim),
  n_genes_pre_filter = n_genes_pre, n_genes_post_filter = n_genes_post,
  n_features_pre_filter = n_feat_pre, n_features_post_filter = n_feat_post,
  common_precision = obs_out$common_precision,
  n_genes_precision_nan_repaired = obs_out$n_genes_precision_repaired,
  n_bins = opts$n_bins,
  min_cells_per_unit = opts$min_cells_per_unit,
  min_units_per_bin = opts$min_units_per_bin,
  min_bin_expr = opts$min_bin_expr,
  keep_induced_false_only = opts$keep_induced_false_only,
  covariate_col = opts$covariate_col,
  master_seed = opts$seed,
  n_sig_gene_pseudotime_bh05 = sig_bh05, n_sig_gene_pseudotime_bh01 = sig_bh01,
  ks_uniformity_p = ks_p,
  start_time = run_start_time, end_time = format(Sys.time()),
  host = run_host
)
summary_path <- file.path(opts$output_dir, sprintf("drimseq_summary_%s.csv", slug))
write.csv(as.data.frame(summary_row), summary_path, row.names = FALSE)
cat("\n[Summary]", summary_path, "\n")
print(as.data.frame(summary_row))

# ---- Sidecar: run_metadata.json ----
metadata <- list(
  script = "260526-drimseq_pseudotime_pseudobulk.R",
  script_version_hash = tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  ),
  R_version = R.version.string,
  package_versions = list(
    anndata = as.character(packageVersion("anndata")),
    DRIMSeq = as.character(packageVersion("DRIMSeq")),
    BiocParallel = as.character(packageVersion("BiocParallel"))
  ),
  master_seed = opts$seed,
  pseudotime_col = opts$pseudotime_col,
  covariate_col = opts$covariate_col,
  lineage = opts$lineage,
  n_bins = opts$n_bins,
  min_cells_per_unit = opts$min_cells_per_unit,
  min_units_per_bin = opts$min_units_per_bin,
  min_bin_expr = opts$min_bin_expr,
  keep_induced_false_only = opts$keep_induced_false_only,
  ncores = opts$ncores,
  start_time = run_start_time,
  end_time = format(Sys.time()),
  host = run_host,
  command = paste(commandArgs(trailingOnly = FALSE), collapse = " "),
  output_gene_csv = gene_csv
)
metadata_path <- file.path(opts$output_dir, sprintf("run_metadata_%s.json", slug))
writeLines(jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), metadata_path)
cat("[Metadata]", metadata_path, "\n")

cat("\n=== Done at", format(Sys.time()), "===\n")
