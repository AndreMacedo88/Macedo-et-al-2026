"""Annotate an h5ad with cell-level metadata from Zenodo metadata.tsv.

Only .obs is modified; .X is read and written unchanged so sparse matrices
remain sparse on disk. Works on both:
  - Cell Ranger-aggregated h5ads (obs column 'sample', obs_name '<BC>-<N>')
  - utrcalling-merged h5ads       (obs column 'sample_id', obs_name '<BC>-1_<N>')
"""
import argparse

import anndata as ad
import pandas as pd
import scipy.sparse as sp


METADATA_COLS = [
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mito",
    "percent.ribo",
    "gRNA",
    "induced",
    "celltype",
]
CATEGORICAL_COLS = {"gRNA", "induced", "celltype"}


def _normalize_library(lib: str) -> str:
    # Convert metadata.tsv's library value to the pipeline's library_id:
    #   "137116"          → "137116"   (Batch 1)
    #   "BATCH3_213664"   → "213664"   (Batch 2 — Zenodo mislabels these)
    #   "206287_rnacount" → "206287"   (Batch 3)
    return lib.removesuffix("_rnacount").removeprefix("BATCH3_")


def _extract_bc_core(cell: str, lib: str) -> str:
    prefix = f"{lib}_"
    trailing = cell[len(prefix):] if cell.startswith(prefix) else cell
    return trailing.rsplit("-", 1)[0]


def _build_metadata_lookup(metadata_path: str) -> pd.DataFrame:
    meta = pd.read_csv(metadata_path, sep="\t")
    libs = meta["library"].astype(str)
    cells = meta["cell"].astype(str)
    meta["_lib_id"] = libs.map(_normalize_library)
    meta["_bc_core"] = [_extract_bc_core(c, l) for c, l in zip(cells, libs)]
    meta = meta.drop_duplicates(subset=["_lib_id", "_bc_core"], keep="first")
    return meta.set_index(["_lib_id", "_bc_core"])


def _extract_keys(adata: ad.AnnData) -> pd.MultiIndex:
    obs = adata.obs
    if "sample_id" in obs.columns:
        lib = obs["sample_id"].astype(str)
    elif "sample" in obs.columns:
        lib = obs["sample"].astype(str)
    else:
        raise ValueError(
            "adata.obs must have 'sample' (Cell Ranger) or 'sample_id' (UTR merged)"
        )
    names = pd.Series(obs.index.astype(str), index=obs.index)
    bc_core = names.str.split("_").str[0].str.rsplit("-", n=1).str[0]
    return pd.MultiIndex.from_arrays([lib.values, bc_core.values])


def annotate(h5ad_path: str, metadata_path: str, output_path: str) -> None:
    print(f"Reading {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)
    print(adata)
    print(f"  X is sparse: {sp.issparse(adata.X)}  ({type(adata.X).__name__})")

    meta = _build_metadata_lookup(metadata_path)
    keys = _extract_keys(adata)
    available = [c for c in METADATA_COLS if c in meta.columns]
    annotations = meta[available].reindex(keys)
    annotations.index = adata.obs.index

    for col in available:
        adata.obs[col] = annotations[col].values
        if col in CATEGORICAL_COLS:
            adata.obs[col] = adata.obs[col].astype("category")

    matched = adata.obs["celltype"].notna().sum() if "celltype" in adata.obs.columns else -1
    print(f"Annotated {matched}/{adata.n_obs} cells with celltype")

    print(f"Writing {output_path}")
    adata.write_h5ad(output_path, compression="gzip")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input h5ad file")
    parser.add_argument("--metadata", required=True, help="Zenodo metadata.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output h5ad file")
    args = parser.parse_args()
    annotate(args.input, args.metadata, args.output)


if __name__ == "__main__":
    main()
