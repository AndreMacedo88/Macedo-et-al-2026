"""Build samples.tsv from the E-MTAB-13148 SDRF.

Fetches the SDRF, keeps only GEX rows (drops AMPLICON libraries — those are
gRNA / CHOOSE-clone-barcode amplicons used for CRISPR demultiplexing, not
full-transcriptome data), groups by library, cross-references
metadata.tsv to determine the Zenodo library name, and writes a samples.tsv
with semicolon-separated ERR and FASTQ URL lists so multi-lane libraries
(Batch 1) round-trip through the pipeline.
"""
import argparse
import re
import urllib.request

import pandas as pd


SDRF_URL = (
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/148/E-MTAB-13148/Files/"
    "E-MTAB-13148.sdrf.txt"
)

STIMULUS_MAP = {
    "induced eCas9, 4OHT treatment": "4OHT_induced",
    "non-induced eCas9, 4OHT treatment": "4OHT_non_induced",
    "no 4OHT treatment": "no_4OHT",
}


def fetch_sdrf(url):
    print(f"Fetching SDRF: {url}")
    with urllib.request.urlopen(url) as r:
        return pd.read_csv(r, sep="\t")


def parse_library_id(assay_name):
    # e.g. sc_b1_137116_137116_S10_L001 → 137116
    #      sc_b1_137122_137122merged_S4_L001 → 137122
    m = re.match(r"sc_b[123]_(\d+)_", assay_name)
    if not m:
        raise ValueError(f"Could not parse library ID from Assay Name: {assay_name}")
    return m.group(1)


def parse_batch(description):
    m = re.search(r"Batch\s+(\d)", description)
    if not m:
        raise ValueError(f"Could not parse batch from Description: {description}")
    return f"batch{m.group(1)}"


def parse_lib_name(description):
    # "GEX, Batch 1, Lib 1" → "GEX_Batch1_Lib1"
    return re.sub(r",\s+", "_", description).replace(" ", "")


def find_metadata_library(lib_id, meta_libs):
    if lib_id in meta_libs:
        return lib_id
    if f"{lib_id}_rnacount" in meta_libs:
        return f"{lib_id}_rnacount"
    if f"BATCH3_{lib_id}" in meta_libs:
        return f"BATCH3_{lib_id}"
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-m", "--metadata", required=True, help="Zenodo metadata.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output samples.tsv")
    parser.add_argument("--sdrf-url", default=SDRF_URL)
    args = parser.parse_args()

    meta_libs = set(
        pd.read_csv(args.metadata, sep="\t", usecols=["library"])["library"]
        .astype(str)
        .unique()
    )
    print(f"Metadata has {len(meta_libs)} unique libraries")

    sdrf = fetch_sdrf(args.sdrf_url)
    print(f"SDRF: {len(sdrf)} rows")

    gex = sdrf[sdrf["Comment[LIBRARY_STRATEGY]"] == "ssRNA-seq"].copy()
    print(f"  GEX rows (ssRNA-seq): {len(gex)}")
    print(f"  Dropped AMPLICON rows: {len(sdrf) - len(gex)}")

    # Keep only cDNA reads (R1 = cell barcode + UMI, R2 = cDNA).
    # Cell Ranger 8.x does not use I1/I2 for standard 10x 3' v3.
    gex["filename"] = gex["Comment[SUBMITTED_FILE_NAME]"].astype(str)
    read_pat = re.compile(r"_(R1|R2)_001\.fastq\.gz$")
    gex = gex[gex["filename"].apply(lambda f: bool(read_pat.search(f)))].copy()
    print(f"  After dropping I1/I2: {len(gex)} rows (R1+R2 only)")

    gex["library_id"] = gex["Assay Name"].astype(str).apply(parse_library_id)
    gex["batch"] = gex["Description"].astype(str).apply(parse_batch)
    gex["library_name"] = gex["Description"].astype(str).apply(parse_lib_name)
    gex["stimulus"] = gex["Characteristics[stimulus]"].map(STIMULUS_MAP)
    gex["err"] = gex["Comment[ENA_RUN]"].astype(str)
    gex["url"] = gex["Comment[FASTQ_URI]"].astype(str)

    unmapped = gex[gex["stimulus"].isna()]
    if len(unmapped) > 0:
        raise ValueError(
            f"Unknown stimulus values: "
            f"{unmapped['Characteristics[stimulus]'].unique().tolist()}"
        )

    rows = []
    warnings = []
    for lib_id, group in gex.groupby("library_id", sort=False):
        group = group.sort_values("filename")
        errs = sorted(group["err"].unique())
        urls = group["url"].tolist()

        metadata_library = find_metadata_library(lib_id, meta_libs)
        if metadata_library is None:
            warnings.append(f"library {lib_id}: no match in metadata.tsv; skipping")
            continue

        rows.append({
            "library_id": lib_id,
            "metadata_library": metadata_library,
            "library_name": group["library_name"].iloc[0],
            "stimulus": group["stimulus"].iloc[0],
            "batch": group["batch"].iloc[0],
            "err_accessions": ";".join(errs),
            "fastq_urls": ";".join(urls),
        })

    for w in warnings:
        print(f"  WARNING: {w}")

    out = pd.DataFrame(rows).sort_values(["batch", "library_id"]).reset_index(drop=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {len(out)} rows to {args.output}")
    print(out[["library_id", "metadata_library", "library_name",
               "stimulus", "batch", "err_accessions"]].to_string(index=False))


if __name__ == "__main__":
    main()
