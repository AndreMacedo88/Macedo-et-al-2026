import argparse
from warnings import simplefilter
import numpy as np
from scipy import sparse
import pandas as pd
import scanpy as sc

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def main():
    parser = argparse.ArgumentParser(
        description="Memory-efficient concatenation of UTR mapping h5ad files."
    )
    parser.add_argument("inputs", nargs="+", help="Input h5ad files")
    parser.add_argument("-o", "--output", required=True, help="Output h5ad file")
    args = parser.parse_args()

    print("---------------------------------")
    print("Reading AnnData objects from file")
    adatas = [sc.read_h5ad(adata) for adata in args.inputs]

    var_set = set()
    var_list = []
    ct = 1
    for adata in adatas:
        print("---------------------------------")
        print(adata)

        adata.obs.index = adata.obs.index + "_" + str(ct)
        adata.obs["batch_number"] = int(ct)

        var = adata.var
        var_index = var.index.to_list()
        var_set.update(var_index)

        var_list.append(var)

        adata.X = sparse.csr_matrix(adata.X)
        ct += 1

    print("---------------------------------")
    print("Creating one single AnnData object and storing")
    adata_concat = sc.concat(adatas=adatas, join="outer")
    print(adata_concat)

    var = pd.concat(var_list, axis=0)
    del var_list
    var = var.drop_duplicates()
    var = var.reindex(adata_concat.var_names.to_list())

    adata_concat.var = var

    print("---------------------------------")
    print("Some Checks")
    print(adata_concat)
    intersection = var_set.intersection(set(adata_concat.var_names))
    pa_n = len(var_set)
    print(f"The new data intersects in {len(intersection)} polyA sites, out of {pa_n}")

    adata_concat.write_h5ad(args.output)


if __name__ == "__main__":
    main()
