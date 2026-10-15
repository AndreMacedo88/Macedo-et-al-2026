"""Write table S1 (GO enrichment behind Fig. 1F) as an Excel file with a caption sheet.

Reads the CSV written by table_s1_go_enrichment.R and writes table_s1_go_enrichment.xlsx next
to it (default folder: $M26_RESULTS_ROOT/tables, or <repo>/results/tables).

Run: conda run -n apa-python python table_s1_xlsx.py [--outdir DIR]
"""

import argparse
import csv
import os
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

CAPTION = (
    "Table S1. Gene Ontology biological processes tested for over-representation among genes "
    "whose poly(A) site usage changes along adult neurogenesis (Fig. 1F)."
)
DETAILS = (
    "Foreground: the 1,746 genes significant (FDR < 5%) under both the asymptotic likelihood-ratio "
    "test and the within-sample permutation test, of which 1,699 have a Gene Ontology "
    "biological-process annotation. Background: the 7,214 tested genes with a converged regression "
    "fit and a Gene Ontology annotation. All 1,436 terms with 50 to 500 annotated background genes "
    "are listed, ordered by P value (one-sided hypergeometric test; clusterProfiler v4.14.0, "
    "org.Mm.eg.db v3.20.0); adjusted P values are Benjamini-Hochberg-adjusted over all 1,436 terms. "
    "The six terms shown in Fig. 1F are marked."
)
# (CSV column, header in the Excel table, description)
COLUMNS: list[tuple[str, str, str]] = [
    ("go_id", "GO ID", "Gene Ontology identifier"),
    ("term", "Biological process", "Gene Ontology term name"),
    ("shown_in_fig1F", "Shown in Fig. 1F", "yes for the six terms drawn in Fig. 1F"),
    ("fg_in_term", "Foreground genes in term", "foreground genes annotated to the term"),
    (
        "fg_annotated",
        "Foreground genes annotated",
        "foreground genes with any biological-process annotation",
    ),
    ("bg_in_term", "Background genes in term", "background genes annotated to the term"),
    (
        "bg_annotated",
        "Background genes annotated",
        "background genes with any biological-process annotation",
    ),
    (
        "fold_enrichment",
        "Fold enrichment",
        "(foreground in term / foreground annotated) / (background in term / background annotated)",
    ),
    ("p_value", "P", "one-sided hypergeometric P value"),
    ("p_adjusted_BH", "Adjusted P", "Benjamini-Hochberg-adjusted P value over all 1,436 terms"),
    ("genes", "Genes", "foreground genes annotated to the term (mouse gene symbols, alphabetical)"),
]
INT_COLS = {"fg_in_term", "fg_annotated", "bg_in_term", "bg_annotated"}
FLOAT_FORMATS = {"fold_enrichment": "0.00", "p_value": "0.00E+00", "p_adjusted_BH": "0.00E+00"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    default_dir = Path(os.environ.get("M26_RESULTS_ROOT", "results")) / "tables"
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", type=Path, default=default_dir)
    outdir = parser.parse_args().outdir
    rows = read_rows(outdir / "table_s1_go_enrichment.csv")
    wb = Workbook()

    info = wb.active
    info.title = "Caption"
    info["A1"] = CAPTION
    info["A1"].font = Font(bold=True)
    info["A2"] = DETAILS
    info["A2"].alignment = Alignment(wrap_text=True, vertical="top")
    info.column_dimensions["A"].width = 30
    info.column_dimensions["B"].width = 90
    info.merge_cells("A2:B2")
    info.row_dimensions[2].height = 120
    info["A4"] = "Column"
    info["B4"] = "Description"
    info["A4"].font = info["B4"].font = Font(bold=True)
    for i, (_, header, desc) in enumerate(COLUMNS, start=5):
        info.cell(row=i, column=1, value=header)
        info.cell(row=i, column=2, value=desc)

    ws = wb.create_sheet("Table S1")
    for j, (_, header, _) in enumerate(COLUMNS, start=1):
        ws.cell(row=1, column=j, value=header).font = Font(bold=True)
    for i, row in enumerate(rows, start=2):
        for j, (col, _, _) in enumerate(COLUMNS, start=1):
            value: str | int | float = row[col]
            if col in INT_COLS:
                value = int(value)
            elif col in FLOAT_FORMATS:
                value = float(value)
            cell = ws.cell(row=i, column=j, value=value)
            if col in FLOAT_FORMATS:
                cell.number_format = FLOAT_FORMATS[col]
    widths = [12, 45, 9, 11, 11, 11, 11, 10, 10, 10, 120]
    for j, w in enumerate(widths, start=1):
        ws.column_dimensions[ws.cell(row=1, column=j).column_letter].width = w
    ws.freeze_panes = "C2"
    xlsx_path = outdir / "table_s1_go_enrichment.xlsx"
    wb.save(xlsx_path)
    print(f"{len(rows)} terms written to {xlsx_path}")


if __name__ == "__main__":
    main()
