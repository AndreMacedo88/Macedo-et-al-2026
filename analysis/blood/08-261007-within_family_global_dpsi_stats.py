#!/usr/bin/env python
"""Summary statistics for the within-family global 3'UTR lengthening (blood, Fig. 4 text).

Purpose
    Archive the text statistics for the within-family proband-vs-sibling global
    Delta-psi result: of the 9 same-batch, age-matched ASD-proband vs unaffected-
    sibling pairs, how many show global lengthening (global_dpsi > 0), the
    two-sided exact Wilcoxon signed-rank P, the two-sided exact sign (binomial)
    test P, and the Spearman correlation between the signed within-pair age
    difference (p_age - s_age) and global_dpsi (age-gap control).

Supports
    Fig. 4 text numbers: 8/9 pairs lengthened, Wilcoxon P = 0.0078,
    sign test P = 0.039, Spearman rho ~ 0.03 (P ~ 0.95) for age gap vs global_dpsi.

Inputs
    $M26_DATA_ROOT/3utr/Blood/2025_analyses/within_family_systematic/per_pair_global_dpsi.csv
    (columns fam, proband, sib, p_age, s_age, samesex, global_dpsi; written by
    07-260615-within_family_systematic_apa.ipynb)

Outputs
    <outdir>/global_dpsi_summary_stats.csv  (one row per statistic: statistic, value, note)
    Default outdir = the input directory (2025_analyses/within_family_systematic/).
    An existing output file is never overwritten unless --overwrite is given.

How to run (conda env apa-python)
    conda run -n apa-python python 08-261007-within_family_global_dpsi_stats.py
    conda run -n apa-python python 08-261007-within_family_global_dpsi_stats.py --outdir /some/scratch/dir

Provenance
    Reconstructed on 2026-10-07 to reproduce the within-family global-Delta-psi text
    statistics (new file global_dpsi_summary_stats.csv); original generating code was
    not archived. The Wilcoxon and sign tests mirror the calls in
    07-260615-within_family_systematic_apa.ipynb (scipy.stats.wilcoxon default,
    scipy.stats.binomtest); the age-gap Spearman is computed on the signed gap
    p_age - s_age.
"""

import argparse
from pathlib import Path

import pandas as pd
from scipy import stats
import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

RES = Path(f"{DATA_ROOT}/3utr/Blood/2025_analyses/within_family_systematic")
IN_CSV = RES / "per_pair_global_dpsi.csv"
OUT_NAME = "global_dpsi_summary_stats.csv"


def summarise(pairs: pd.DataFrame) -> pd.DataFrame:
    """Compute the global-Delta-psi summary statistics for the sibling pairs."""
    g = pairs["global_dpsi"].to_numpy(dtype=float)
    n = len(g)
    n_pos = int((g > 0).sum())
    wil = stats.wilcoxon(g, alternative="two-sided", method="exact")
    sign = stats.binomtest(n_pos, n, p=0.5, alternative="two-sided")
    age_gap = (pairs["p_age"] - pairs["s_age"]).to_numpy(dtype=float)
    rho, rho_p = stats.spearmanr(age_gap, g)
    rows = [
        ("n_pairs", n, "same-batch ASD-proband vs unaffected-sibling pairs"),
        ("n_global_dpsi_gt0", n_pos, "pairs with global_dpsi > 0 (proband longer)"),
        ("median_global_dpsi", float(pd.Series(g).median()), "median over pairs"),
        ("mean_global_dpsi", float(g.mean()), "mean over pairs"),
        (
            "wilcoxon_signed_rank_stat",
            float(wil.statistic),
            "scipy W (smaller rank sum)",
        ),
        ("wilcoxon_signed_rank_p", float(wil.pvalue), "two-sided exact, H0 median = 0"),
        (
            "sign_test_p",
            float(sign.pvalue),
            f"two-sided exact binomial, {n_pos}/{n} > 0, p0 = 0.5",
        ),
        (
            "spearman_rho_agegap_vs_global_dpsi",
            float(rho),
            "age gap = p_age - s_age (signed, years)",
        ),
        (
            "spearman_p_agegap_vs_global_dpsi",
            float(rho_p),
            "two-sided, scipy.stats.spearmanr",
        ),
    ]
    return pd.DataFrame(rows, columns=["statistic", "value", "note"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", type=Path, default=IN_CSV, help="per_pair_global_dpsi.csv")
    ap.add_argument("--outdir", type=Path, default=RES, help="output directory")
    ap.add_argument("--overwrite", action="store_true", help="allow replacing an existing output")
    args = ap.parse_args()

    pairs = pd.read_csv(args.input)
    tab = summarise(pairs)
    print(tab.to_string(index=False))

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / OUT_NAME
    if out.exists() and not args.overwrite:
        print(f"{out} exists - not overwritten (use --overwrite)")
        return
    tab.to_csv(out, index=False)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
