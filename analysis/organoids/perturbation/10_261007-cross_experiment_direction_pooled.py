"""Pooled direction agreement between treatment and developmental 3'UTR shifts (Fig. 2G text).

Purpose
    For the gene-lineage pairs significant in both experiments (cross-experiment overlap:
    sig_treat & sig_pseudo in overlap_gene_table_<lin>.csv; 25 INH + 12 EXC = 37 pairs),
    test whether the covariate-adjusted treatment Delta-psi and the developmental
    (pseudotime) delta_psi share the same sign, per lineage and pooled:
    same-sign count, two-sided exact binomial test vs 0.5, Pearson and Spearman correlation.

Supports
    Fig. 2G text statistics (pooled: 15/37 same sign, binomial P = 0.324;
    Pearson r = -0.181, P = 0.284; per lineage INH 10/25, EXC 5/12).

Provenance
    Written on 2026-10-07; reproduces the pooled direction statistics reported in the paper
    and stores them in overlap_direction_pooled.csv. Per-lineage rows reproduce the ADJ rows of
    cross_experiment_overlap/overlap_direction_agreement_adjusted.csv.

Inputs (KB = $M26_DATA_ROOT/3utr/organoids/knoblich/utrcalling)
    KB/dirmult_treatment_results/cross_experiment_overlap/overlap_gene_table_<lin>.csv
    KB/dirmult_treatment_results/<lin>/main_only/delta_psi_treatment_adjusted_<lin>.csv
        (from 04_261007-make_adjusted_dpsi.py; column delta_psi_treatment_adjusted)
    KB/composite_ranking/composite_ranking_optimized2_<lin>.csv (column delta_psi, matched
        on gene_name == overlap symbol)

Output (default KB/dirmult_treatment_results/cross_experiment_overlap/)
    overlap_direction_pooled.csv: one row each for inhibitory, excitatory, pooled with
    n_pairs, same_sign, sign_concordance, binom_p_two_sided, pearson_r, pearson_p,
    spearman_rho, spearman_p

Run (conda env apa-python)
    conda run -n apa-python python 10_261007-cross_experiment_direction_pooled.py [--outdir DIR]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

KB = Path(f"{DATA_ROOT}/3utr/organoids/knoblich/utrcalling")
TREAT_DIR = KB / "dirmult_treatment_results"
OVERLAP_DIR = TREAT_DIR / "cross_experiment_overlap"
LINEAGES = ["inhibitory", "excitatory"]


def overlap_pairs(lin: str) -> pd.DataFrame:
    ov = pd.read_csv(
        OVERLAP_DIR / f"overlap_gene_table_{lin}.csv", keep_default_na=False
    )
    ov = ov[
        ov["sig_treat"].astype(str).eq("True") & ov["sig_pseudo"].astype(str).eq("True")
    ]
    adj = pd.read_csv(
        TREAT_DIR / lin / "main_only" / f"delta_psi_treatment_adjusted_{lin}.csv",
        keep_default_na=False,
    )
    comp = pd.read_csv(
        KB / "composite_ranking" / f"composite_ranking_optimized2_{lin}.csv",
        keep_default_na=False,
    )
    dev = comp.loc[comp["gene_name"].isin(ov["symbol"]), ["gene_name", "delta_psi"]]
    assert dev["gene_name"].is_unique, f"[{lin}] ambiguous composite gene_name match"

    pairs = (
        ov[["symbol"]]
        .merge(
            adj[["gene_id", "delta_psi_treatment_adjusted"]],
            left_on="symbol",
            right_on="gene_id",
        )
        .merge(
            dev.rename(columns={"delta_psi": "delta_psi_dev"}),
            left_on="symbol",
            right_on="gene_name",
        )
    )
    assert len(pairs) == len(ov), (
        f"[{lin}] {len(ov) - len(pairs)} overlap genes not matched"
    )
    pairs["lineage"] = lin
    return pairs[["lineage", "symbol", "delta_psi_treatment_adjusted", "delta_psi_dev"]]


def direction_stats(df: pd.DataFrame, scope: str) -> dict[str, float | int | str]:
    x = df["delta_psi_treatment_adjusted"].astype(float).to_numpy()
    y = df["delta_psi_dev"].astype(float).to_numpy()
    n = len(df)
    same = int((np.sign(x) == np.sign(y)).sum())
    pear = stats.pearsonr(x, y)
    spear = stats.spearmanr(x, y)
    return {
        "scope": scope,
        "n_pairs": n,
        "same_sign": same,
        "sign_concordance": same / n,
        "binom_p_two_sided": stats.binomtest(same, n, 0.5).pvalue,
        "pearson_r": pear.statistic,
        "pearson_p": pear.pvalue,
        "spearman_rho": spear.statistic,
        "spearman_p": spear.pvalue,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--outdir", type=Path, default=OVERLAP_DIR)
    args = parser.parse_args()

    pairs = pd.concat([overlap_pairs(lin) for lin in LINEAGES], ignore_index=True)
    rows = [direction_stats(pairs[pairs["lineage"] == lin], lin) for lin in LINEAGES]
    rows.append(direction_stats(pairs, "pooled"))
    res = pd.DataFrame(rows)

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / "overlap_direction_pooled.csv"
    res.to_csv(out, index=False)
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"-> {out}")


if __name__ == "__main__":
    main()
