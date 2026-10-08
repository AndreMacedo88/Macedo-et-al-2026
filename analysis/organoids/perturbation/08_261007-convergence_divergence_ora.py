"""Convergence / divergence of treatment ORA terms between EXC and INH lineages.

Purpose
    Merge the per-lineage treatment ORA tables (one tier) and classify every tested term:
        Convergent = p.adjust < 0.05 in both lineages
        EXC-only   = p.adjust < 0.05 in excitatory only
        INH-only   = p.adjust < 0.05 in inhibitory only
        NS         = otherwise (including terms tested in one lineage only and not significant)
    Collections (in this order): Reactome, GO_BP, Hallmark, KEGG.

Supports
    Fig. 2H (convergent / lineage-specific perturbation pathways) and fig. S3H.

Provenance
    Ported on 2026-10-07 from an earlier R script that wrote
    dirmult_treatment_results/convergence_divergence_ora.csv. The port keeps the R merge
    semantics (dplyr::full_join on ID: EXC rows in EXC order, then INH-only-tested
    terms in INH order; Description taken from the EXC table only, so it is NA for terms tested
    only in INH) and R write.csv formatting (quoted strings, 15 significant digits, NA).

Inputs (KB = $M26_DATA_ROOT/3utr/organoids/knoblich/utrcalling/dirmult_treatment_results)
    KB/<lin>/enrichment_<tier>/ora_{Reactome,GO_BP,Hallmark,KEGG}.csv   (lin = excitatory,
    inhibitory; written by 07_260515-enrichment_discovery.ipynb for the discovery tier, which
    the figures use; the robust tier is the same analysis on the bootstrap-robust subset)

Outputs (default directory KB)
    --tier robust    -> convergence_divergence_ora.csv
    --tier discovery -> convergence_divergence_ora_discovery.csv   (default; used for the figures)
    columns: ID, Description, padj_exc, k_exc, padj_inh, k_inh, collection, class

Run (conda env apa-python)
    conda run -n apa-python python 08_261007-convergence_divergence_ora.py [--tier discovery|robust]
        [--outdir DIR]

Verification (2026-10-07): --tier robust output is byte-identical to the archived
    convergence_divergence_ora.csv (46 Convergent, 109 EXC-only, 6 INH-only).
"""

import argparse
import math
from pathlib import Path

import pandas as pd
import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

KB = Path(
    f"{DATA_ROOT}/3utr/organoids/knoblich/utrcalling/dirmult_treatment_results"
)
COLLECTIONS = ["Reactome", "GO_BP", "Hallmark", "KEGG"]
PADJ_CUT = 0.05
OUT_NAMES = {
    "robust": "convergence_divergence_ora.csv",
    "discovery": "convergence_divergence_ora_discovery.csv",
}


def read_ora(tier: str, lin: str, coll: str) -> pd.DataFrame:
    path = KB / lin / f"enrichment_{tier}" / f"ora_{coll}.csv"
    return pd.read_csv(path, float_precision="round_trip")


def merge_collection(tier: str, coll: str) -> pd.DataFrame:
    """dplyr::full_join(exc[ID, Description, padj_exc, k_exc], inh[ID, padj_inh, k_inh], by = 'ID')."""
    exc = read_ora(tier, "excitatory", coll).rename(
        columns={"p.adjust": "padj_exc", "Count": "k_exc"}
    )
    inh = read_ora(tier, "inhibitory", coll).rename(
        columns={"p.adjust": "padj_inh", "Count": "k_inh"}
    )
    exc = exc[["ID", "Description", "padj_exc", "k_exc"]]
    inh = inh[["ID", "padj_inh", "k_inh"]]
    assert exc["ID"].is_unique and inh["ID"].is_unique, f"duplicated IDs in {coll}"
    left = exc.merge(inh, on="ID", how="left")
    right_only = inh[~inh["ID"].isin(exc["ID"])]
    merged = pd.concat([left, right_only], ignore_index=True)
    merged["collection"] = coll

    sig_exc = merged["padj_exc"].notna() & (merged["padj_exc"] < PADJ_CUT)
    sig_inh = merged["padj_inh"].notna() & (merged["padj_inh"] < PADJ_CUT)
    merged["class"] = "NS"
    merged.loc[sig_exc & ~sig_inh, "class"] = "EXC-only"
    merged.loc[sig_inh & ~sig_exc, "class"] = "INH-only"
    merged.loc[sig_exc & sig_inh, "class"] = "Convergent"
    return merged


def r_double(x: float) -> str:
    """Format a double like R write.csv (15 significant digits, scipen = 0)."""
    if pd.isna(x):
        return "NA"
    if x == 0:
        return "0"
    sig = f"{x:.15g}"
    mantissa = sig.split("e")[0].lstrip("-").replace(".", "").lstrip("0")
    n_digits = max(len(mantissa.rstrip("0")), 1)
    exponent = math.floor(math.log10(abs(float(sig))))
    sci = f"{x:.{n_digits - 1}e}"
    decimals = max(0, n_digits - 1 - exponent)
    fixed = f"{x:.{decimals}f}"
    return fixed if len(fixed) <= len(sci) else sci


def r_int(x: float) -> str:
    return "NA" if pd.isna(x) else str(int(x))


def r_str(x: object) -> str:
    if pd.isna(x):
        return "NA"
    return '"' + str(x).replace('"', '""') + '"'


def write_r_csv(df: pd.DataFrame, path: Path) -> None:
    fmt = {
        "ID": r_str,
        "Description": r_str,
        "padj_exc": r_double,
        "k_exc": r_int,
        "padj_inh": r_double,
        "k_inh": r_int,
        "collection": r_str,
        "class": r_str,
    }
    lines = [",".join(f'"{c}"' for c in df.columns)]
    for row in df.itertuples(index=False):
        lines.append(",".join(fmt[c](v) for c, v in zip(df.columns, row)))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tier", choices=["robust", "discovery"], default="discovery")
    parser.add_argument("--outdir", type=Path, default=KB)
    args = parser.parse_args()

    cd = pd.concat(
        [merge_collection(args.tier, coll) for coll in COLLECTIONS], ignore_index=True
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / OUT_NAMES[args.tier]
    write_r_csv(cd, out)

    counts = cd["class"].value_counts()
    print(f"[{args.tier}] {len(cd)} terms -> {out}")
    print(
        "  "
        + " | ".join(
            f"{cl}: {counts.get(cl, 0)}"
            for cl in ["Convergent", "EXC-only", "INH-only"]
        )
    )
    print(cd.groupby(["collection", "class"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
