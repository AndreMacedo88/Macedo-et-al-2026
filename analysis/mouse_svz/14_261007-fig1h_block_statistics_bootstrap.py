"""Block statistics and gene-bootstrap CIs for the bin-by-bin similarity matrices of Fig. 1H.

Purpose:
    Fig. 1H shows Spearman similarity between 20 pseudotime bins (qNSC1 -> NB, adult SVZ)
    computed from the transcriptome (log CPM per gene) and from the 3'UTRome (expression-
    weighted relative 3'UTR length W per gene). This script quantifies, per layer, the mean
    min-max-normalised similarity across two cell-type boundaries (TAP x NB and
    qNSC2 x aNSC1) and within the NB block, the difference 3'UTRome minus transcriptome, and
    percentile 95% CIs from B = 500 gene-bootstrap replicates (numpy default_rng seed 0).

Supports: Fig. 1H numbers quoted in the Results text (e.g. TAP x NB 0.776 vs 0.631,
    diff -0.144 [-0.208, -0.066]; within-NB 0.920 vs 0.773, diff -0.147 [-0.193, -0.106]).

Algorithm (same as the original analysis):
    1. Cells of the core lineage (qNSC1, qNSC2, aNSC1, aNSC2, TAP, NB) with finite DPT
       pseudotime are split into 20 equal-count pseudotime bins (quantile edges).
    2. Per bin, raw 3'UTR counts are pseudobulked. For each gene: W = sum over its 3'UTR
       isoforms of (isoform fraction x within-gene rank of utr_size scaled to [0, 1]);
       E = log1p(gene counts / bin library size x 1e6). Bins with < 20 gene reads are NaN.
    3. Genes covered in all bins; top 1000 most variable genes across bins per layer.
    4. Spearman correlation between bins (genes ranked within each bin), min-max normalised
       using the off-diagonal values of each layer's matrix.
    5. Block statistics on the bin majority cell-type labels; gene bootstrap resamples the
       top-1000 genes of each layer independently (W first, then E, per replicate).

Inputs (same as 13_260701-correlations_over_pseudotime.ipynb):
    utrcalling/UTR_count_table_genelevel_bw12_above10cts_apa_named_adata_harmony_lineage.h5ad
        (cell types, dpt_pseudotime)
    utrcalling/UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad
        (raw 3'UTR isoform counts, var 'name' and 'utr_size')

Output:
    utrcalling/correlation_analyses_proportions/fig1h_block_statistics_bootstrap.csv
    (tidy: one row per statistic x layer, with point estimate and percentile 95% CI)

How to run:
    conda run -n apa-python python 14_261007-fig1h_block_statistics_bootstrap.py [--outdir DIR] [--overwrite]

Ported on 2026-10-07 from earlier analysis scripts (2026-09-03) that printed the Fig. 1H
block statistics but wrote no file.
"""

import argparse
import warnings
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import rankdata
import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

warnings.filterwarnings("ignore")

UTR_DIR = Path(f"{DATA_ROOT}/3utr/wildtypes_analysis/2026_analyses/utrcalling")
LINEAGE_H5AD = UTR_DIR / "UTR_count_table_genelevel_bw12_above10cts_apa_named_adata_harmony_lineage.h5ad"
RAW_H5AD = UTR_DIR / "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad"
DEFAULT_OUTDIR = UTR_DIR / "correlation_analyses_proportions"
OUT_NAME = "fig1h_block_statistics_bootstrap.csv"

CORE = ["qNSC1", "qNSC2", "aNSC1", "aNSC2", "TAP", "NB"]
N_BINS = 20
MIN_READS = 20
TOP_N = 1000
N_BOOT = 500
SEED = 0
# Bin majority labels as used in the original analysis (checked against the data below).
EXPECTED_MAJORITY = ["qNSC1"] * 4 + ["qNSC2"] + ["aNSC1"] * 2 + ["aNSC2"] * 2 + ["TAP"] * 2 + ["NB"] * 9
STATISTICS = {
    "tap_x_nb": ("TAP", "NB", "mean between-block similarity, TAP bins x NB bins"),
    "qnsc2_x_ansc1": ("qNSC2", "aNSC1", "mean between-block similarity, qNSC2 bins x aNSC1 bins"),
    "nb_within": ("NB", "NB", "mean within-block similarity over all NB bin pairs (i<j)"),
}


def build_bin_profiles() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return W (3'UTRome) and E (transcriptome) gene x bin matrices and bin majority labels."""
    lineage = ad.read_h5ad(LINEAGE_H5AD)
    raw = ad.read_h5ad(RAW_H5AD)[lineage.obs_names].copy()
    celltype = lineage.obs["celltype"].astype(str).values
    pseudotime = lineage.obs["dpt_pseudotime"].astype(float).values
    X = raw.X.tocsr() if issparse(raw.X) else raw.X
    site_gene = raw.var["name"].astype(str).values
    site_size = raw.var["utr_size"].astype(float).values

    genes = list(dict.fromkeys(site_gene))
    sites_by_gene: dict[str, list[int]] = defaultdict(list)
    for i, g in enumerate(site_gene):
        sites_by_gene[g].append(i)
    gene_sites = {g: np.array(v) for g, v in sites_by_gene.items()}
    # within-gene rank of 3'UTR length, scaled to [0, 1] (single-isoform genes -> 0)
    rel_rank = {
        g: ((rankdata(site_size[i], method="average") - 1) / (len(i) - 1) if len(i) > 1 else np.zeros(len(i)))
        for g, i in gene_sites.items()
    }

    core = np.isin(celltype, CORE) & np.isfinite(pseudotime)
    pt_core = pseudotime[core]
    X_core = X[core]
    edges = np.quantile(pt_core, np.linspace(0, 1, N_BINS + 1))
    edges[0] -= 1e-9
    bin_of = np.clip(np.digitize(pt_core, edges) - 1, 0, N_BINS - 1)
    majority = np.array(
        [pd.Series(celltype[core][bin_of == b]).value_counts().idxmax() for b in range(N_BINS)]
    )

    bin_sums = np.vstack([np.asarray(X_core[bin_of == b].sum(axis=0)).ravel() for b in range(N_BINS)])
    lib_size = bin_sums.sum(axis=1)
    W = np.full((len(genes), N_BINS), np.nan)
    E = np.full((len(genes), N_BINS), np.nan)
    for gi, g in enumerate(genes):
        sub = bin_sums[:, gene_sites[g]]
        gene_sum = sub.sum(axis=1)
        ok = gene_sum >= MIN_READS
        if ok.any():
            W[gi, ok] = (sub[ok] / gene_sum[ok, None]) @ rel_rank[g]
            E[gi, ok] = np.log1p(gene_sum[ok] / lib_size[ok] * 1e6)
    return W, E, majority


def top_variable(P: np.ndarray, covered: np.ndarray, n: int) -> np.ndarray:
    """Indices of the n covered genes with the highest variance across bins."""
    v = np.full(P.shape[0], np.nan)
    v[covered] = np.nanvar(P[covered], axis=1)
    order = np.argsort(v)[::-1]
    return order[np.isfinite(v[order])][:n]


def spearman_bins(P: np.ndarray) -> np.ndarray:
    """Bin x bin Spearman correlation of a gene x bin matrix."""
    R = np.apply_along_axis(rankdata, 0, P)
    return np.corrcoef(R.T)


def block_stats(C: np.ndarray, blocks: dict[str, np.ndarray]) -> dict[str, float]:
    """Min-max normalise C by its off-diagonal range, then average over each block pair."""
    v = C[np.triu_indices_from(C, k=1)]
    M = (C - v.min()) / (v.max() - v.min())
    out = {}
    for key, (a, b, _) in STATISTICS.items():
        sub = M[np.ix_(blocks[a], blocks[b])]
        out[key] = float(sub[np.triu_indices(len(blocks[a]), k=1)].mean() if a == b else sub.mean())
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--overwrite", action="store_true", help="replace an existing output CSV")
    args = parser.parse_args()
    out_path = args.outdir / OUT_NAME
    if out_path.exists() and not args.overwrite:
        raise SystemExit(f"{out_path} exists; pass --overwrite to replace it")

    W, E, majority = build_bin_profiles()
    if list(majority) != EXPECTED_MAJORITY:
        raise SystemExit(f"bin majority labels changed: {list(majority)}")
    blocks = {c: np.where(majority == c)[0] for c in CORE}

    covered = np.isfinite(W).all(axis=1) & np.isfinite(E).all(axis=1)
    top_w = top_variable(W, covered, TOP_N)
    top_e = top_variable(E, covered, TOP_N)
    obs_w = block_stats(spearman_bins(W[top_w]), blocks)
    obs_e = block_stats(spearman_bins(E[top_e]), blocks)

    rng = np.random.default_rng(SEED)
    boot: dict[str, dict[str, list[float]]] = {k: {"w": [], "e": []} for k in STATISTICS}
    for _ in range(N_BOOT):
        idx_w = rng.choice(len(top_w), len(top_w), replace=True)
        idx_e = rng.choice(len(top_e), len(top_e), replace=True)
        bw = block_stats(spearman_bins(W[top_w[idx_w]]), blocks)
        be = block_stats(spearman_bins(E[top_e[idx_e]]), blocks)
        for k in STATISTICS:
            boot[k]["w"].append(bw[k])
            boot[k]["e"].append(be[k])

    rows = []
    for k, (a, b, desc) in STATISTICS.items():
        reps = {
            "transcriptome": np.array(boot[k]["e"]),
            "utrome": np.array(boot[k]["w"]),
            "utrome_minus_transcriptome": np.array(boot[k]["w"]) - np.array(boot[k]["e"]),
        }
        point = {
            "transcriptome": obs_e[k],
            "utrome": obs_w[k],
            "utrome_minus_transcriptome": obs_w[k] - obs_e[k],
        }
        for layer, est in point.items():
            lo, hi = np.percentile(reps[layer], [2.5, 97.5])
            rows.append(
                dict(statistic=k, block_a=a, block_b=b, layer=layer, estimate=est,
                     ci95_low=lo, ci95_high=hi, n_bins_a=len(blocks[a]), n_bins_b=len(blocks[b]),
                     n_top_genes=TOP_N, n_boot=N_BOOT, seed=SEED, description=desc)
            )
    table = pd.DataFrame(rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path, index=False)
    with pd.option_context("display.width", 200):
        print(table[["statistic", "layer", "estimate", "ci95_low", "ci95_high"]].round(3).to_string(index=False))
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
