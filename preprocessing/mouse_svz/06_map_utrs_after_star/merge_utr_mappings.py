import os
import click
from warnings import simplefilter
import numpy as np
from scipy import sparse
import pandas as pd
import scanpy as sc

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


@click.command()
@click.argument("inputs", nargs=-1)
@click.option("-o", "--output", help="Path to write the results", required=True)
def main(inputs, output):
    print("---------------------------------")
    print("Reading AnnData objects from file")
    adatas = [sc.read_h5ad(adata) for adata in inputs]

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

    adata_concat.write_h5ad(output)


if __name__ == "__main__":
    main()
