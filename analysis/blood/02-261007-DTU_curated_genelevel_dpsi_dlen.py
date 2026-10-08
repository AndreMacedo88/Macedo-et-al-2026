#!/usr/bin/env python
"""Per-gene ASD-vs-CTRL Delta-psi and delta-UTR-length on the curated blood cohort.

Purpose
    Effect-size companions to the gene-level DRIMSeq DTU
    (01-261007-DTU_curated_genelevel_unadj.R) for the curated 71-sample cohort
    (obs status_split in {Train, Test_original}; 40 ASD / 31 CTRL):
      * gene_dpsi_curated.csv  - Delta-psi: per sample and gene,
            W = sum_k x_k * r_k / sum_k x_k,   r_k = (rank(site position) - 1) / (K - 1)
        (x_k = raw site counts, K = number of sites of the gene, position = utr_size);
        dpsi = mean W(ASD) - mean W(CTRL), samples with zero counts for the gene excluded.
      * gene_dlen_curated.csv  - the same in bp: W_len = proportion-weighted mean site
        position, W_len = sum_k p_k * L_k / sum_k p_k with p_k = x_k / sum_k x_k
        (float32 proportions, i.e. exactly the values of the
        UTR_count_table_..._apa_named_percentages.h5ad layer); dlen = mean W_len(ASD) -
        mean W_len(CTRL), zero-count samples excluded. Positive = longer / more distal in ASD.
    All 8,338 genes are reported (2 genes are NaN: no counts in one diagnosis group).

Supports
    Fig. 4B (volcano x-axis and direction of the 81 DTU genes) and fig. S4C; also read by
    figures/fig4/fig5c_convergence.ipynb and analysis/blood/07-260615-within_family_systematic_apa.ipynb.

Inputs
    $M26_DATA_ROOT/3utr/Blood/all_till_20240724/utrcalling/
        UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad

Outputs (default outdir = $M26_DATA_ROOT/3utr/Blood/2025_analyses/DTU_curated_sexadj)
    gene_dpsi_curated.csv (gene_id, dpsi), gene_dlen_curated.csv (gene_id, dlen).
    Existing files are never overwritten unless --overwrite is given.

How to run (conda env apa-python; < 1 min)
    conda run -n apa-python python 02-261007-DTU_curated_genelevel_dpsi_dlen.py --outdir /path/to/scratch

Provenance
    Reconstructed on 2026-10-07 to reproduce DTU_curated_sexadj/gene_dpsi_curated.csv
    (2026-06-15) and gene_dlen_curated.csv (2026-06-14); original generating code was
    not archived. The W computation follows 07-260615-within_family_systematic_apa.ipynb.
    Verification (2026-10-07): both CSVs are byte-identical to the archived files (cmp).
    Bitwise identity depends on the exact float path used here: float32 proportions with
    float32 per-gene sums for dlen, and selecting group columns by sample label (not a
    boolean mask) before .mean(axis=1); other summation orders differ by <= 1e-12.
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

H5 = Path(
    f"{DATA_ROOT}/3utr/Blood/all_till_20240724/utrcalling/"
    "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad"
)
OUTDIR = Path(f"{DATA_ROOT}/3utr/Blood/2025_analyses/DTU_curated_sexadj")
CURATED = ["Train", "Test_original"]


def rank_weights(gene: np.ndarray, length: np.ndarray) -> np.ndarray:
    """r_k = (rank of site position within gene - 1) / (K - 1); 0 for single-site genes."""
    r = np.zeros(len(gene))
    for _, idx in pd.DataFrame({"g": gene}).groupby("g").indices.items():
        k = len(idx)
        r[idx] = 0.0 if k < 2 else (rankdata(length[idx], method="average") - 1) / (k - 1)
    return r


def group_delta(w: pd.DataFrame, asd: np.ndarray, ctrl: np.ndarray) -> pd.Series:
    """mean over ASD columns minus mean over CTRL columns (NaN = excluded sample)."""
    return w.loc[:, asd].mean(axis=1) - w.loc[:, ctrl].mean(axis=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--h5ad", type=Path, default=H5)
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    a = ad.read_h5ad(args.h5ad)
    x = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)  # float32 counts
    gene = a.var["utr_name"].astype(str).to_numpy()
    length = pd.to_numeric(a.var["utr_size"]).to_numpy().astype(float)
    samples = a.obs_names.to_numpy()
    cur = a.obs["status_split"].astype(str).isin(CURATED).to_numpy()
    dx = a.obs["Diagnosis"].astype(str).to_numpy()
    asd = samples[cur & (dx == "ASD")]
    ctrl = samples[cur & (dx == "CTRL")]
    print(f"curated n = {cur.sum()} (ASD {len(asd)} / CTRL {len(ctrl)}); genes = {len(set(gene))}")

    # Delta-psi from raw counts (float64 numerator, float32 integer denominator)
    r = rank_weights(gene, length)
    num = pd.DataFrame((x * r[None, :]).T, index=gene).groupby(level=0).sum()
    den = pd.DataFrame(x.T, index=gene).groupby(level=0).sum()
    w = num / den.replace(0, np.nan)
    w.columns = samples
    dpsi = group_delta(w, asd, ctrl).rename("dpsi")

    # delta-length in bp from float32 within-gene proportions
    den_site = den.loc[gene].to_numpy().T  # per-site gene total, samples x sites (float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        prop = np.where(den_site > 0, x / den_site, np.float32(0)).astype(np.float32)
    num_len = pd.DataFrame((prop * length[None, :]).T, index=gene).groupby(level=0).sum()
    den_prop = pd.DataFrame(prop.T, index=gene).groupby(level=0).sum()  # float32 sums
    w_len = num_len / den_prop.replace(0, np.nan)
    w_len.columns = samples
    dlen = group_delta(w_len, asd, ctrl).rename("dlen")

    args.outdir.mkdir(parents=True, exist_ok=True)
    for series, name in [
        (dpsi, "gene_dpsi_curated.csv"),
        (dlen, "gene_dlen_curated.csv"),
    ]:
        out = args.outdir / name
        print(
            f"{name}: {len(series)} genes, NaN {int(series.isna().sum())}, >0 {int((series > 0).sum())}"
        )
        if out.exists() and not args.overwrite:
            print(f"  {out} exists - not overwritten (use --overwrite)")
            continue
        series.rename_axis("gene_id").to_csv(out)
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
