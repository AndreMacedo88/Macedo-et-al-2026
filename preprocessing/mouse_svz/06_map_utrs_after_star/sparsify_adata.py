import click
import scanpy as sc
from scipy import sparse


@click.command()
@click.argument("input")
@click.argument("output")
def main(input, output):

    adata = sc.read_h5ad(input)
    adata.X = sparse.csr_matrix(adata.X)
    adata.write_h5ad(filename=output, compression="gzip")


if __name__ == "__main__":
    main()
