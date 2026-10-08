"""Aggregate per-sample Cell Ranger filtered matrices into a single AnnData.

Merges the matrices, rewrites per-sample barcode suffixes (-N) Cell Ranger
aggr-style, and uses Ensembl gene IDs as var_names. Cell-level metadata is
not added here; the `annotate_cells` rule handles that downstream so there's
one canonical join point.
"""

import argparse

import anndata as ad
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Per-sample filtered_feature_bc_matrix.h5 files; order defines the -N suffix.",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample/library IDs, one per input, same order.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output .h5ad")
    args = parser.parse_args()

    if len(args.inputs) != len(args.samples):
        raise ValueError(
            f"--samples has {len(args.samples)} entries but got {len(args.inputs)} inputs"
        )

    adatas = []
    for i, (h5, sample) in enumerate(zip(args.inputs, args.samples), start=1):
        print(f"[{i}/{len(args.inputs)}] Reading {sample} from {h5}")
        adata = sc.read_10x_h5(h5)
        adata.var["gene_symbol"] = adata.var_names.astype(str)
        adata.var_names = adata.var["gene_ids"].astype(str)
        adata.var_names.name = "gene_id"
        del adata.var["gene_ids"]
        adata.obs_names = adata.obs_names.str.replace(r"-\d+$", f"-{i}", regex=True)
        adata.obs["sample"] = sample
        adatas.append(adata)

    print("Concatenating...")
    adata_concat = ad.concat(adatas, join="outer", merge="same")
    adata_concat.obs["sample"] = adata_concat.obs["sample"].astype("category")
    print(adata_concat)

    print(f"Writing to {args.output}")
    adata_concat.write_h5ad(args.output)


if __name__ == "__main__":
    main()
