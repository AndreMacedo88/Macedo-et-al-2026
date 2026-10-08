"""Model-matched, covariate-adjusted treatment Delta-psi per gene (Fig. 2F x-axis).

Purpose
    For every gene, compute the 'induced' coefficient of a weighted least-squares
    regression of the per-unit 3'UTR length index W_u on the DRIMSeq treatment-LRT
    design  W_u ~ bs(pt_u, df=5) + super_batch + induced  (weights = n_cells), i.e. the
    treatment effect on 3'UTR length holding pseudotime and super_batch fixed. Also
    recomputes the unadjusted (cell-count-weighted marginal) Delta-psi on the same units.

Supports
    Fig. 2F (treatment volcano, x-axis = delta_psi_treatment_adjusted).

Provenance
    Rewritten on 2026-10-07; reproduces the stored
    dirmult_treatment_results/<lin>/main_only/delta_psi_treatment_adjusted_<lin>.csv.

Definitions
    Units: the pseudobulk units of the main DRIMSeq run
        (02_drimseq_treatment_main_pseudobulk.R --require_targeting_guide TRUE): cells with
        non-NA dpt_pseudotime_<lin>; induced cells whose gRNA is 'Control2' or NA removed;
        pseudotime quantile-binned into 8 bins over the retained cells (R quantile type 7,
        cut(include.lowest=TRUE)); unit = sample_id x bin; units with >= 10 cells kept;
        pt_mid = mean pseudotime of all retained cells in the bin; super_batch = first 3
        characters of sample_id. Matches main_only/pseudobulk_units.csv (checked at run time).
    Sites: all var features of a gene (grouped by var.gene_name, var order = proximal ->
        distal); r_k = (k-1)/(K-1); pi_uk = Y_uk / sum_j Y_uj; W_u = sum_k r_k pi_uk,
        computed on units with > 0 counts for the gene.
    Genes kept: >= 11 expressed units (= number of design columns + 2) with at least one
        induced and one uninduced unit.
    Raw Delta-psi: n_cells-weighted mean W (induced) - n_cells-weighted mean W (uninduced).
    Adjusted Delta-psi: WLS via numpy lstsq on sqrt(n_cells)-scaled patsy design
        'bs(pt_mid, df=5) + C(super_batch) + induced' (spline knots from all units of the
        lineage, all super_batch levels kept even if absent for a gene, hence pinv covariance);
        sigma2 = resid'resid / (n_units - n_columns); SE from pinv(A'A); two-sided t-test p
        (adj_wls_p); BH across genes with scipy.stats.false_discovery_control
        (adj_wls_padj). Output sorted by gene_id (gene symbol).

Inputs (KB = $M26_DATA_ROOT/3utr/organoids/knoblich/utrcalling)
    KB/
        merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad
        dirmult_treatment_results/<lin>/main_only/pseudobulk_units.csv  (consistency check only)

Outputs (default: KB/dirmult_treatment_results/<lin>/main_only/)
    delta_psi_treatment_adjusted_<lin>.csv with columns gene_id,
    delta_psi_treatment_adjusted, delta_psi_treatment_raw_recomputed, adj_wls_p, adj_wls_padj

Run (conda env apa-python; < 1 min, ~3 GB RAM)
    conda run -n apa-python python 04_261007-make_adjusted_dpsi.py [--outdir DIR]
        [--lineages inhibitory excitatory]
    With --outdir, files are written to DIR/<lin>/.

Verification (2026-10-07): units identical to pseudobulk_units.csv; output CSVs are
    byte-identical to the archived files for both lineages (numpy/scipy/patsy in env apa-python).
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import patsy
import scipy.sparse as sp
from scipy import stats
import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

KB = Path(f"{DATA_ROOT}/3utr/organoids/knoblich/utrcalling")
H5AD = KB / "merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad"
TREAT_DIR = KB / "dirmult_treatment_results"
CONTROL_GUIDE = "Control2"
N_BINS = 8
MIN_CELLS_PER_UNIT = 10
SPLINE_DF = 5
DESIGN = f"bs(pt_mid, df={SPLINE_DF}) + C(super_batch) + induced_i"


def quantile_bin(pt: np.ndarray, n_bins: int) -> np.ndarray:
    """Replicates R quantile(type 7) + cut(include.lowest=TRUE, labels=FALSE)."""
    qs = np.unique(np.quantile(pt, np.linspace(0, 1, n_bins + 1)))
    bins = np.searchsorted(qs, pt, side="left")
    bins[bins == 0] = 1
    return bins


def build_units(adata: ad.AnnData, lin: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Pseudobulk (sample x pseudotime-bin) units as in 02_drimseq_treatment_main_pseudobulk.R."""
    obs = adata.obs
    pt = obs[f"dpt_pseudotime_{lin}"].to_numpy(dtype=float)
    induced = obs["induced"].astype(str).isin(["True", "TRUE", "true", "1"]).to_numpy()
    grna = obs["gRNA"].astype(object).to_numpy()
    drop_ctrl = induced & (pd.isna(grna) | (grna == CONTROL_GUIDE))
    keep = ~np.isnan(pt) & ~drop_ctrl

    pt_k = pt[keep]
    sample_k = obs["sample_id"].astype(str).to_numpy()[keep]
    bins = quantile_bin(pt_k, N_BINS)
    bin_mid = {b: pt_k[bins == b].mean() for b in np.unique(bins)}
    groups = np.array([f"{s}__bin{b:02d}" for s, b in zip(sample_k, bins)])
    unit_keys, inv = np.unique(groups, return_inverse=True)

    agg = sp.csr_matrix(
        (np.ones(len(inv)), (inv, np.arange(len(inv)))),
        shape=(len(unit_keys), len(inv)),
    )
    counts = (agg @ adata.X[keep]).toarray()
    induced_by_sample = (
        pd.Series(induced[keep], index=sample_k).groupby(level=0).first()
    )

    meta = pd.DataFrame({"sample_id": unit_keys})
    meta["sample_id_root"] = meta["sample_id"].str.replace(r"__bin\d+$", "", regex=True)
    meta["bin_id"] = meta["sample_id"].str.extract(r"__bin(\d+)$")[0].astype(int)
    meta["pt_mid"] = meta["bin_id"].map(bin_mid)
    meta["induced"] = meta["sample_id_root"].map(induced_by_sample).astype(bool)
    meta["super_batch"] = meta["sample_id_root"].str[:3]
    meta["n_cells"] = np.bincount(inv)
    keep_units = (meta["n_cells"] >= MIN_CELLS_PER_UNIT).to_numpy()
    print(
        f"[{lin}] cells retained {keep.sum()} | units kept {keep_units.sum()} / {len(meta)}"
    )
    return meta[keep_units].reset_index(drop=True), counts[keep_units]


def check_units(meta: pd.DataFrame, lin: str) -> None:
    manifest = TREAT_DIR / lin / "main_only" / "pseudobulk_units.csv"
    if not manifest.exists():
        print(f"[{lin}] no pseudobulk_units.csv to check against")
        return
    ref = pd.read_csv(manifest, dtype={"sample_id_root": str, "super_batch": str})
    m = ref.merge(meta, on="sample_id", suffixes=("_ref", ""))
    ok = (
        len(m) == len(ref) == len(meta)
        and (m["n_cells_ref"] == m["n_cells"]).all()
        and np.allclose(m["pt_mid_ref"], m["pt_mid"], rtol=0, atol=1e-12)
        and (m["super_batch_ref"] == m["super_batch"]).all()
    )
    print(f"[{lin}] units identical to DRIMSeq manifest: {ok}")


def adjusted_dpsi(
    meta: pd.DataFrame, counts: np.ndarray, var: pd.DataFrame
) -> pd.DataFrame:
    design = meta.assign(induced_i=meta["induced"].astype(int))
    x_full = patsy.dmatrix(DESIGN, design, return_type="dataframe")
    i_ind = list(x_full.columns).index("induced_i")
    x_full = x_full.to_numpy()
    n_cols = x_full.shape[1]
    min_units = n_cols + 2
    n_cells = meta["n_cells"].to_numpy(dtype=float)
    is_ind = meta["induced"].to_numpy()

    rows = []
    for gene, idx in var.groupby("gene_name", sort=False).indices.items():
        n_sites = len(idx)
        if n_sites < 2:
            continue
        y = counts[:, idx].astype(float)
        tot = y.sum(axis=1)
        ku = tot > 0
        it = is_ind[ku]
        if ku.sum() < min_units or not it.any() or it.all():
            continue
        rank = np.arange(n_sites) / (n_sites - 1)
        w = (y[ku] / tot[ku, None]) @ rank
        n = n_cells[ku]

        raw = (w[it] * n[it]).sum() / n[it].sum() - (w[~it] * n[~it]).sum() / n[
            ~it
        ].sum()

        sw = np.sqrt(n)
        a = x_full[ku] * sw[:, None]
        yb = w * sw
        beta = np.linalg.lstsq(a, yb, rcond=None)[0]
        df_resid = len(n) - n_cols
        resid = yb - a @ beta
        sigma2 = resid @ resid / df_resid
        se = np.sqrt(sigma2 * np.linalg.pinv(a.T @ a)[i_ind, i_ind])
        p = 2 * stats.t.sf(abs(beta[i_ind] / se), df_resid)
        rows.append((gene, beta[i_ind], raw, p))

    out = pd.DataFrame(
        rows,
        columns=[
            "gene_id",
            "delta_psi_treatment_adjusted",
            "delta_psi_treatment_raw_recomputed",
            "adj_wls_p",
        ],
    )
    out["adj_wls_padj"] = stats.false_discovery_control(
        out["adj_wls_p"].to_numpy(), method="bh"
    )
    return out.sort_values("gene_id").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Write to OUTDIR/<lin>/ instead of the canonical main_only/ directories",
    )
    parser.add_argument("--lineages", nargs="+", default=["inhibitory", "excitatory"])
    args = parser.parse_args()

    adata = ad.read_h5ad(H5AD)
    adata.X = adata.X.tocsr()
    var = adata.var[["gene_name"]].astype(str).reset_index(drop=True)

    for lin in args.lineages:
        meta, counts = build_units(adata, lin)
        check_units(meta, lin)
        res = adjusted_dpsi(meta, counts, var)
        outdir = (
            TREAT_DIR / lin / "main_only" if args.outdir is None else args.outdir / lin
        )
        outdir.mkdir(parents=True, exist_ok=True)
        out = outdir / f"delta_psi_treatment_adjusted_{lin}.csv"
        res.to_csv(out, index=False)
        print(
            f"[{lin}] {len(res)} genes | adj_wls_padj < 0.05: {(res['adj_wls_padj'] < 0.05).sum()} -> {out}"
        )


if __name__ == "__main__":
    main()
