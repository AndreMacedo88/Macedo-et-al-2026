#!/usr/bin/env python3
"""
Delta Psi (Δψ) — Rank-Weighted Proportion Shift for 3'UTR Isoform Usage.

Adapts the PSI (Percent Spliced In) framework to APA: for each gene with K UTR
sites ordered by utr_size, assigns normalized ranks r_k = (k-1)/(K-1), fits
the same GLM model as the DirMult R script (per-site weighted logistic GLM +
softmax), and computes a signed scalar in [-1, 1] summarizing the direction
and magnitude of UTR length change along pseudotime.

Effect size metrics per gene:
  - delta_psi:     W(t_late) - W(t_early)  — overall endpoint change
  - delta_psi_max: max(W) - min(W)         — max trajectory range
  - nmi:           1 - |delta_psi|/|delta_psi_max| — non-monotonicity index in [0, 1]
  - t_at_max/min:  where W(t) peaks/troughs along pseudotime

Usage:
  python 260309-calculate_delta_psi.py <input.h5ad> [options]
"""

import argparse
import os
import sys
import warnings
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import anndata
import numpy as np
import pandas as pd
import patsy
import scipy.sparse as sp
import statsmodels.api as sm
from scipy.stats import rankdata
from tqdm import tqdm

warnings.filterwarnings("ignore", category=sm.tools.sm_exceptions.ConvergenceWarning)
warnings.filterwarnings(
    "ignore", category=sm.tools.sm_exceptions.PerfectSeparationWarning
)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# GLM fitting and prediction
# ---------------------------------------------------------------------------


def fit_gene_glms(counts, design_matrix, n_sites):
    """Fit K independent weighted logistic GLMs (mirrors R rough_dirmnom_ll).

    Parameters
    ----------
    counts : ndarray (n_cells, K)
        Raw UMI counts for expressing cells.
    design_matrix : ndarray (n_cells, p)
        Full design matrix (batch + spline terms).
    n_sites : int
        Number of UTR sites (K).

    Returns
    -------
    beta : ndarray (p, K)
        Coefficient matrix.
    """
    row_sums = counts.sum(axis=1)
    fracs = counts / row_sums[:, None]
    p = design_matrix.shape[1]
    beta = np.zeros((p, n_sites))

    for j in range(n_sites):
        y = fracs[:, j]
        try:
            glm = sm.GLM(
                y,
                design_matrix,
                family=sm.families.Binomial(link=sm.families.links.Logit()),
                freq_weights=row_sums,
            )
            result = glm.fit(maxiter=50, disp=False)
            beta[:, j] = result.params
        except Exception:
            # Fallback: use log-proportion as intercept, zeros elsewhere
            mean_frac = np.clip(y.mean(), 1e-8, 1 - 1e-8)
            beta[0, j] = np.log(mean_frac / (1 - mean_frac))
    return beta


def predict_proportions(pred_matrix, beta):
    """Softmax prediction: π̂ = softmax(X @ β).

    Parameters
    ----------
    pred_matrix : ndarray (n_grid, p)
    beta : ndarray (p, K)

    Returns
    -------
    pi_hat : ndarray (n_grid, K)
    """
    eta = pred_matrix @ beta
    # Numerical stability: subtract row max before exp
    eta -= eta.max(axis=1, keepdims=True)
    exp_eta = np.exp(eta)
    pi_hat = exp_eta / exp_eta.sum(axis=1, keepdims=True)
    return pi_hat


# ---------------------------------------------------------------------------
# Per-gene worker
# ---------------------------------------------------------------------------


def process_gene(args):
    """Process a single gene: fit GLM, predict trajectory, compute Δψ.

    Designed to be called via multiprocessing.Pool.map.
    """
    (
        gene_id,
        gene_name,
        col_indices,
        utr_sizes,
        site_names,
        X_csc,
        design_full,
        design_nobatch,
        pred_rows_per_batch,
        batch_weights,
        batch_ids,
        cell_batches,
        t_grid,
        n_grid,
        min_counts,
        min_cells,
        has_batch,
        ranks,
    ) = args

    n_sites = len(col_indices)
    result_gene = {
        "gene_id": gene_id,
        "gene_name": gene_name,
        "n_sites": n_sites,
        "n_cells": np.nan,
        "total_umis": np.nan,
        "delta_psi": np.nan,
        "delta_psi_max": np.nan,
        "nmi": np.nan,
        "t_at_min": np.nan,
        "t_at_max": np.nan,
        "tvd": np.nan,
        "max_abs_delta_pi": np.nan,
    }
    result_sites = []

    try:
        # Extract sparse slice -> dense (small: n_cells x K)
        gene_counts = X_csc[:, col_indices].toarray()
        row_sums = gene_counts.sum(axis=1)
        expressing = row_sums > 0
        n_expr = expressing.sum()
        total_umis = row_sums[expressing].sum()

        result_gene["n_cells"] = int(n_expr)
        result_gene["total_umis"] = float(total_umis)

        if n_expr < min_cells or total_umis < min_counts:
            # Return NaN result + empty site results
            for j in range(n_sites):
                result_sites.append(
                    {
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "site_name": site_names[j],
                        "utr_size": utr_sizes[j],
                        "rank": ranks[j],
                        "pi_early": np.nan,
                        "pi_late": np.nan,
                        "delta_pi": np.nan,
                    }
                )
            return result_gene, result_sites

        counts_expr = gene_counts[expressing]
        batches_expr = cell_batches[expressing]
        active_batches = np.unique(batches_expr)

        # Pick design matrix: with or without batch term
        if has_batch and len(active_batches) > 1:
            dm = design_full[expressing]
            use_batch = True
        else:
            dm = design_nobatch[expressing]
            use_batch = False

        # Fit GLMs
        beta = fit_gene_glms(counts_expr, dm, n_sites)

        # Predict proportions at grid points, averaged across batches
        pi_trajectory = np.zeros((n_grid, n_sites))
        total_weight = 0.0

        for b in active_batches:
            if use_batch:
                pred_mat = pred_rows_per_batch[b]
            else:
                # No batch term: use the no-batch prediction matrix
                pred_mat = pred_rows_per_batch["nobatch"]
            pi_b = predict_proportions(pred_mat, beta)
            w = batch_weights.get(b, 1.0)
            pi_trajectory += pi_b * w
            total_weight += w

        if total_weight > 0:
            pi_trajectory /= total_weight

        # Compute W(t) = sum_k r_k * pi_k(t)
        W = pi_trajectory @ ranks

        # Effect sizes
        delta_psi = W[-1] - W[0]
        W_min_idx = np.argmin(W)
        W_max_idx = np.argmax(W)
        delta_psi_max = W[W_max_idx] - W[W_min_idx]
        # Sign delta_psi_max: positive if max comes after min1
        if W_max_idx < W_min_idx:
            delta_psi_max = -delta_psi_max
        t_at_min = t_grid[W_min_idx]
        t_at_max = t_grid[W_max_idx]

        # Per-site endpoint proportions
        pi_early = pi_trajectory[0]
        pi_late = pi_trajectory[-1]
        delta_pi = pi_late - pi_early
        tvd = 0.5 * np.abs(delta_pi).sum()
        max_abs_delta_pi = np.abs(delta_pi).max()

        # NMI = 1 - |Δψ|/|Δψ_max| in [0, 1]: 0 = monotone, →1 = transient/bell-shaped
        if abs(delta_psi_max) < 1e-12:
            nmi = 0.0
        else:
            nmi = float(np.clip(1.0 - abs(delta_psi) / abs(delta_psi_max), 0.0, 1.0))

        result_gene.update(
            {
                "delta_psi": float(np.clip(delta_psi, -1, 1)),
                "delta_psi_max": float(np.clip(delta_psi_max, -1, 1)),
                "nmi": nmi,
                "t_at_min": float(t_at_min),
                "t_at_max": float(t_at_max),
                "tvd": float(tvd),
                "max_abs_delta_pi": float(max_abs_delta_pi),
            }
        )

        for j in range(n_sites):
            result_sites.append(
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "site_name": site_names[j],
                    "utr_size": utr_sizes[j],
                    "rank": ranks[j],
                    "pi_early": float(pi_early[j]),
                    "pi_late": float(pi_late[j]),
                    "delta_pi": float(delta_pi[j]),
                }
            )

    except Exception as e:
        # GLM convergence or numerical failure — return NaN
        for j in range(n_sites):
            result_sites.append(
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "site_name": site_names[j],
                    "utr_size": utr_sizes[j],
                    "rank": ranks[j] if j < len(ranks) else np.nan,
                    "pi_early": np.nan,
                    "pi_late": np.nan,
                    "delta_pi": np.nan,
                }
            )

    return result_gene, result_sites


# ---------------------------------------------------------------------------
# Design matrix construction
# ---------------------------------------------------------------------------


def build_design_matrices(
    pseudotime, batch, spline_df, t_grid, batch_ids, batch_cell_counts, has_batch
):
    """Build design matrices for fitting and prediction.

    Returns
    -------
    design_full : ndarray (n_cells, p_full)
        Design matrix with batch + spline terms.
    design_nobatch : ndarray (n_cells, p_nobatch)
        Design matrix with spline terms only.
    pred_rows_per_batch : dict
        {batch_id: ndarray (n_grid, p)} for each batch, plus "nobatch" key.
    """
    df = pd.DataFrame(
        {
            "pseudotime": pseudotime,
            "batch": pd.Categorical(batch),
        }
    )

    # Full model: batch + spline
    if has_batch:
        formula_full = f"C(batch) + bs(pseudotime, df={spline_df})"
        design_full = np.asarray(
            patsy.dmatrix(formula_full, data=df, return_type="matrix")
        )
        di_full = patsy.dmatrix(formula_full, data=df).design_info
    else:
        design_full = None
        di_full = None

    # No-batch model: spline only
    formula_nobatch = f"bs(pseudotime, df={spline_df})"
    design_nobatch = np.asarray(
        patsy.dmatrix(formula_nobatch, data=df, return_type="matrix")
    )
    di_nobatch = patsy.dmatrix(formula_nobatch, data=df).design_info

    # Build prediction matrices at grid points
    pred_rows = {}

    if has_batch:
        for b in batch_ids:
            df_pred = pd.DataFrame(
                {
                    "pseudotime": t_grid,
                    "batch": pd.Categorical(
                        [b] * len(t_grid), categories=df["batch"].cat.categories
                    ),
                }
            )
            pred_rows[b] = np.asarray(
                patsy.build_design_matrices([di_full], df_pred)[0]
            )

    # No-batch prediction
    df_pred_nb = pd.DataFrame({"pseudotime": t_grid})
    pred_rows["nobatch"] = np.asarray(
        patsy.build_design_matrices([di_nobatch], df_pred_nb)[0]
    )

    return design_full, design_nobatch, pred_rows


# ---------------------------------------------------------------------------
# Gene index building
# ---------------------------------------------------------------------------


def build_gene_index(var_df, gene_col):
    """Group features by gene, sort by utr_size, compute ranks.

    Returns list of dicts with keys: gene_id, gene_name, col_indices,
    utr_sizes, site_names, ranks.
    """
    genes = defaultdict(list)
    for idx, (feat_name, row) in enumerate(var_df.iterrows()):
        gid = row[gene_col]
        genes[gid].append((idx, feat_name, row))

    gene_list = []
    for gid, entries in genes.items():
        if len(entries) < 2:
            continue
        # Sort by utr_size
        entries.sort(key=lambda x: x[2]["utr_size"])
        col_indices = [e[0] for e in entries]
        utr_sizes = [float(e[2]["utr_size"]) for e in entries]
        site_names = [e[1] for e in entries]
        gene_name = entries[0][2].get("gene_name", "")
        if pd.isna(gene_name):
            gene_name = ""

        K = len(entries)
        if K == 1:
            ranks = np.array([0.0])
        else:
            # Handle tied utr_sizes with midrank
            raw_ranks = rankdata(utr_sizes, method="average")
            ranks = (raw_ranks - 1) / (K - 1)

        gene_list.append(
            {
                "gene_id": gid,
                "gene_name": gene_name,
                "col_indices": col_indices,
                "utr_sizes": utr_sizes,
                "site_names": site_names,
                "ranks": ranks,
            }
        )
    return gene_list


# ---------------------------------------------------------------------------
# Discrete two-group contrast (e.g. neuroblast -> neuron)
# ---------------------------------------------------------------------------


def run_discrete_mode(adata, args, output_dir):
    """Δψ between two DISCRETE cell groups (e.g. NB vs Neurons).

    Rationale: when the two states are disconnected on the manifold (no usable
    continuous pseudotime between them), a continuous spline is indefensible.
    Here Δψ = W(late) - W(early), with W = sum_k r_k * pi_k and r_k the same
    rank-normalised UTR position used in pseudotime mode. pi_k is the EMPIRICAL
    UTR-site usage proportion of a group, pooled within each batch and then
    averaged across batches with EQUAL weight (breaks the batch/group confound
    when the late group is concentrated in a few batches). No GLM, no p-value:
    Δψ is reported as a bounded effect size only.
    """
    if not args.group_early or not args.group_late:
        sys.exit("discrete mode requires --group_early and --group_late")
    early_labels = [s.strip() for s in args.group_early.split(",") if s.strip()]
    late_labels = [s.strip() for s in args.group_late.split(",") if s.strip()]

    ct = adata.obs[args.group_col].astype(str).values
    batch = adata.obs[args.batch_col].astype(str).values
    in_groups = np.isin(ct, early_labels + late_labels)
    if args.restrict_batches.strip():
        keep_batches = [s.strip() for s in args.restrict_batches.split(",") if s.strip()]
        in_batches = np.isin(batch, keep_batches)
    else:
        keep_batches = sorted(set(batch))
        in_batches = np.ones(len(batch), dtype=bool)

    keep = in_groups & in_batches
    adata = adata[keep].copy()
    ct = ct[keep]
    batch = batch[keep]
    is_early = np.isin(ct, early_labels)
    is_late = np.isin(ct, late_labels)

    print(
        f"  discrete contrast: early={early_labels} (n={int(is_early.sum())}) "
        f"vs late={late_labels} (n={int(is_late.sum())}); batches kept={keep_batches}"
    )
    grp = np.where(is_late, "late", np.where(is_early, "early", "other"))
    print(pd.crosstab(pd.Series(grp, name="group"), pd.Series(batch, name="batch")).to_string())

    X = adata.X
    X_csc = X.tocsc() if sp.issparse(X) else sp.csc_matrix(X)
    gene_list = build_gene_index(adata.var, args.gene_col)
    print(
        f"  {len(gene_list)} multi-site genes; gate = >={args.min_cells} cells "
        f"AND >={args.min_counts} UMIs in EACH group"
    )
    batches_present = [b for b in keep_batches if (batch == b).any()]

    def group_pi(cnt, grp_mask):
        """Equal-batch-weighted empirical proportions for one group."""
        pis = []
        ncell = 0
        umis = 0.0
        for b in batches_present:
            m = grp_mask & (batch == b)
            if not m.any():
                continue
            sub = cnt[m]
            subrow = sub.sum(axis=1)
            ncell += int((subrow > 0).sum())
            umis += float(subrow.sum())
            tot = np.asarray(sub.sum(axis=0)).ravel()
            s = tot.sum()
            if s > 0:
                pis.append(tot / s)
        if not pis:
            return None, ncell, umis
        return np.mean(np.vstack(pis), axis=0), ncell, umis

    gene_results = []
    site_results = []
    for g in tqdm(gene_list, desc="Genes(discrete)"):
        col_idx = g["col_indices"]
        ranks = np.asarray(g["ranks"])
        K = len(col_idx)
        cnt = X_csc[:, col_idx].toarray()
        pi_e, nce, ue = group_pi(cnt, is_early)
        pi_l, ncl, ul = group_pi(cnt, is_late)
        rec = {
            "gene_id": g["gene_id"],
            "gene_name": g["gene_name"],
            "n_sites": K,
            "n_cells_early": nce,
            "n_cells_late": ncl,
            "umis_early": ue,
            "umis_late": ul,
            "n_cells": nce + ncl,
            "total_umis": ue + ul,
            "W_early": np.nan,
            "W_late": np.nan,
            "delta_psi": np.nan,
        }
        ok = (
            pi_e is not None
            and pi_l is not None
            and nce >= args.min_cells
            and ncl >= args.min_cells
            and ue >= args.min_counts
            and ul >= args.min_counts
        )
        if ok:
            We = float(ranks @ pi_e)
            Wl = float(ranks @ pi_l)
            rec["W_early"] = We
            rec["W_late"] = Wl
            rec["delta_psi"] = float(np.clip(Wl - We, -1, 1))
        gene_results.append(rec)
        for j in range(K):
            site_results.append(
                {
                    "gene_id": g["gene_id"],
                    "gene_name": g["gene_name"],
                    "site_name": g["site_names"][j],
                    "utr_size": g["utr_sizes"][j],
                    "rank": float(ranks[j]),
                    "pi_early": float(pi_e[j]) if ok else np.nan,
                    "pi_late": float(pi_l[j]) if ok else np.nan,
                    "delta_pi": float(pi_l[j] - pi_e[j]) if ok else np.nan,
                }
            )

    tag = args.output_tag or (
        "-".join(early_labels) + "_to_" + "-".join(late_labels)
    ).lower()
    df_gene = pd.DataFrame(gene_results).sort_values(
        "delta_psi", key=lambda s: s.abs(), ascending=False, na_position="last"
    )
    df_site = pd.DataFrame(site_results).sort_values(["gene_id", "utr_size"])
    gene_path = os.path.join(output_dir, f"delta_psi_{tag}_per_gene.csv")
    site_path = os.path.join(output_dir, f"delta_psi_{tag}_per_site.csv")
    df_gene.to_csv(gene_path, index=False)
    df_site.to_csv(site_path, index=False)

    valid = df_gene["delta_psi"].notna()
    print(f"\n[discrete:{tag}] {int(valid.sum())}/{len(df_gene)} genes pass gates")
    if valid.any():
        d = df_gene.loc[valid, "delta_psi"]
        print(f"  Δψ range [{d.min():.4f}, {d.max():.4f}], median|Δψ| = {d.abs().median():.4f}")
        print(f"  lengthening (+): {int((d > 0).sum())}   shortening (-): {int((d < 0).sum())}")
    print(f"  Saved: {gene_path}")
    print(f"         {site_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Compute Δψ (delta psi) effect sizes for 3'UTR isoform usage along pseudotime."
    )
    parser.add_argument("input_path", help="Path to h5ad with raw counts")
    parser.add_argument(
        "--spline_df",
        "-d",
        type=int,
        default=5,
        help="B-spline degrees of freedom (default: 5)",
    )
    parser.add_argument(
        "--q_early",
        type=float,
        default=0.1,
        help="Early pseudotime quantile (default: 0.1)",
    )
    parser.add_argument(
        "--q_late",
        type=float,
        default=0.9,
        help="Late pseudotime quantile (default: 0.9)",
    )
    parser.add_argument(
        "--n_grid", type=int, default=50, help="Pseudotime grid points (default: 50)"
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        default=100,
        help="Min total UMIs per gene (default: 100)",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        default=30,
        help="Min expressing cells per gene (default: 30)",
    )
    parser.add_argument(
        "--n_cores",
        type=int,
        default=None,
        help="Cores for parallelization (default: cpu_count - 1)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Output directory (default: DirMult/ next to input)",
    )
    parser.add_argument(
        "--batch_weight",
        type=str,
        default="cells",
        choices=["cells", "equal"],
        help="Batch averaging: 'cells' (default) or 'equal'",
    )
    parser.add_argument(
        "--pseudotime_col",
        type=str,
        default="dpt_pseudotime",
        help="obs column for pseudotime (default: dpt_pseudotime)",
    )
    parser.add_argument(
        "--batch_col",
        type=str,
        default="batch_number",
        help="obs column for batch (default: batch_number)",
    )
    parser.add_argument(
        "--gene_col",
        type=str,
        default="utr_name",
        help="var column for gene grouping (default: utr_name)",
    )
    # --- discrete two-group contrast (e.g. neuroblast vs neuron) ---
    parser.add_argument(
        "--mode",
        type=str,
        default="pseudotime",
        choices=["pseudotime", "discrete"],
        help="pseudotime (default; continuous spline) or discrete (two-group contrast)",
    )
    parser.add_argument(
        "--group_col",
        type=str,
        default="celltype",
        help="[discrete] obs column with group labels (default: celltype)",
    )
    parser.add_argument(
        "--group_early",
        type=str,
        default=None,
        help="[discrete] early-state label(s), comma-separated (e.g. 'NB')",
    )
    parser.add_argument(
        "--group_late",
        type=str,
        default=None,
        help="[discrete] late-state label(s), comma-separated (e.g. 'Neurons')",
    )
    parser.add_argument(
        "--restrict_batches",
        type=str,
        default="",
        help="[discrete] comma-separated batches to keep (e.g. '3,4'); empty = all",
    )
    parser.add_argument(
        "--output_tag",
        type=str,
        default=None,
        help="[discrete] filename tag, e.g. 'nb_to_neuron' -> delta_psi_nb_to_neuron_per_gene.csv",
    )
    args = parser.parse_args()

    n_cores = args.n_cores or max(1, cpu_count() - 1)
    output_dir = args.output_dir or os.path.join(
        os.path.dirname(args.input_path), "DirMult"
    )
    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load & validate
    # ------------------------------------------------------------------
    print(f"Loading {args.input_path}")
    adata = anndata.read_h5ad(args.input_path)

    if args.mode == "discrete":
        run_discrete_mode(adata, args, output_dir)
        return

    # Filter NaN pseudotime
    pt = adata.obs[args.pseudotime_col].values.astype(float)
    valid = ~np.isnan(pt)
    if valid.sum() < valid.size:
        print(f"  Dropping {(~valid).sum()} cells with NaN pseudotime")
        adata = adata[valid]
        pt = pt[valid]

    # Validate counts are non-negative integer-like
    X = adata.X
    if sp.issparse(X):
        sample_vals = X.data[: min(10000, len(X.data))]
    else:
        sample_vals = X.ravel()[:10000]
    if np.any(sample_vals < 0):
        print(
            "WARNING: Negative values detected — data may be log-normalized!",
            file=sys.stderr,
        )
    if not np.allclose(sample_vals, np.round(sample_vals), atol=0.01):
        print(
            "WARNING: Non-integer values detected — data may be normalized!",
            file=sys.stderr,
        )

    # Ensure CSC for efficient column slicing
    if sp.issparse(X):
        X_csc = X.tocsc()
    else:
        X_csc = sp.csc_matrix(X)

    batch = adata.obs[args.batch_col].values
    batch_ids = np.unique(batch)
    has_batch = len(batch_ids) > 1

    print(
        f"  {adata.n_obs} cells × {adata.n_vars} features, "
        f"{len(batch_ids)} batch(es)"
    )

    # ------------------------------------------------------------------
    # 2. Pseudotime quantiles and grid
    # ------------------------------------------------------------------
    t_early = float(np.quantile(pt, args.q_early))
    t_late = float(np.quantile(pt, args.q_late))
    t_grid = np.linspace(t_early, t_late, args.n_grid)
    print(
        f"  Pseudotime grid: {t_early:.4f} (q{args.q_early}) → "
        f"{t_late:.4f} (q{args.q_late}), {args.n_grid} points"
    )

    # ------------------------------------------------------------------
    # 3. Build design matrices
    # ------------------------------------------------------------------
    print("Building design matrices...")
    batch_cell_counts = {b: (batch == b).sum() for b in batch_ids}
    if args.batch_weight == "equal":
        batch_weights = {b: 1.0 for b in batch_ids}
    else:
        batch_weights = {b: float(c) for b, c in batch_cell_counts.items()}

    design_full, design_nobatch, pred_rows = build_design_matrices(
        pt, batch, args.spline_df, t_grid, batch_ids, batch_cell_counts, has_batch
    )

    # ------------------------------------------------------------------
    # 4. Build gene index
    # ------------------------------------------------------------------
    print("Building gene index...")
    gene_list = build_gene_index(adata.var, args.gene_col)
    print(f"  {len(gene_list)} multi-site genes")

    # ------------------------------------------------------------------
    # 5. Process genes in parallel
    # ------------------------------------------------------------------
    print(f"Processing genes with {n_cores} cores...")

    # Build argument tuples for worker function
    worker_args = []
    for g in gene_list:
        worker_args.append(
            (
                g["gene_id"],
                g["gene_name"],
                g["col_indices"],
                g["utr_sizes"],
                g["site_names"],
                X_csc,
                design_full,
                design_nobatch,
                pred_rows,
                batch_weights,
                batch_ids,
                batch,
                t_grid,
                args.n_grid,
                args.min_counts,
                args.min_cells,
                has_batch,
                g["ranks"],
            )
        )

    gene_results = []
    site_results = []

    if n_cores == 1:
        for wa in tqdm(worker_args, desc="Genes"):
            rg, rs = process_gene(wa)
            gene_results.append(rg)
            site_results.extend(rs)
    else:
        # Use imap_unordered for progress tracking
        with Pool(n_cores) as pool:
            for rg, rs in tqdm(
                pool.imap_unordered(process_gene, worker_args, chunksize=20),
                total=len(worker_args),
                desc="Genes",
            ):
                gene_results.append(rg)
                site_results.extend(rs)

    # ------------------------------------------------------------------
    # 6. Compile & save
    # ------------------------------------------------------------------
    df_gene = pd.DataFrame(gene_results)
    df_site = pd.DataFrame(site_results)

    # Sort gene results by absolute delta_psi
    df_gene = df_gene.sort_values(
        "delta_psi", key=lambda s: s.abs(), ascending=False, na_position="last"
    )
    df_site = df_site.sort_values(["gene_id", "utr_size"])

    tag = f"df{args.spline_df}_q{args.q_early}-{args.q_late}"
    gene_path = os.path.join(output_dir, f"delta_psi_per_gene_{tag}.csv")
    site_path = os.path.join(output_dir, f"delta_psi_per_site_{tag}.csv")

    df_gene.to_csv(gene_path, index=False)
    df_site.to_csv(site_path, index=False)

    # Summary stats
    valid_genes = df_gene["delta_psi"].notna()
    n_valid = valid_genes.sum()
    print(f"\nResults: {n_valid}/{len(df_gene)} genes with valid Δψ")
    if n_valid > 0:
        dpsi = df_gene.loc[valid_genes, "delta_psi"]
        print(f"  Δψ range: [{dpsi.min():.4f}, {dpsi.max():.4f}]")
        print(f"  Median |Δψ|: {dpsi.abs().median():.4f}")
        print(f"  Shortening (Δψ < 0): {(dpsi < 0).sum()} genes")
        print(f"  Lengthening (Δψ > 0): {(dpsi > 0).sum()} genes")

    print(f"\nSaved:")
    print(f"  Gene-level: {gene_path}")
    print(f"  Site-level:  {site_path}")


if __name__ == "__main__":
    main()
