
import os 
import anndata as ad 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import scanpy as sc
from sklearn.model_selection import train_test_split 

from scipy.sparse import csr_matrix

## VIASH START
par = {
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        'resources_dir': 'src/process_data/',
        'util_dir': 'src/utils',
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

from helper_data import sum_by
from util import download_and_uncompress_zip


def add_metadata(adata, name: str):
    adata.uns['dataset_summary'] = 'Processed experimental single-cell gene expression datasets used in BEELINE.'
    adata.uns['dataset_description'] = 'Processed experimental single-cell gene expression datasets used in BEELINE.'
    adata.uns['data_reference'] = "@dataset{aditya_pratapa_2020_3701939,\nauthor       = {Aditya Pratapa and Amogh Jalihal and Jeffrey Law and Aditya Bharadwaj and T M Murali},\n title        = {Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data },\nmonth        = mar,\nyear         = 2020,\npublisher    = {Zenodo},\ndoi          = {10.5281/zenodo.3701939},\nurl          = {https://doi.org/10.5281/zenodo.3701939},\n}"
    adata.uns['data_url'] = 'https://zenodo.org/records/3701939'
    adata.uns['dataset_id'] = f'beeline_{name}'
    adata.uns['dataset_name'] = f'BEELINE_{name}'
    adata.uns['dataset_organism'] = 'human' if name.startswith('h') else 'mouse'
    adata.uns['normalization_id'] = 'X_norm'
    return adata


if __name__ == '__main__':

    # - get the data
    download_and_uncompress_zip(
        "https://zenodo.org/records/3701939/files/BEELINE-data.zip?download=1",
        "resources/datasets_raw/beeline"
    )

    # Convert CSV files to H5AD
    for name in ["hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L"]:
        filepath = f"resources/datasets_raw/beeline/BEELINE-data/inputs/scRNA-Seq/{name}/ExpressionData.csv"

        df = pd.read_csv(filepath, index_col=0).transpose()
        adata = ad.AnnData(X=df.values, obs=df.index.to_frame(index=False), var=pd.DataFrame(index=df.columns))
        adata.var.index.name = "gene_name"
        adata.layers['X_norm'] = adata.X.copy()

        # Make all obs/var column names strings and clean up a bit
        adata.obs.columns = pd.Index(adata.obs.columns).map(str).str.replace("/", "_")
        adata.var.columns = pd.Index(adata.var.columns).map(str).str.replace("/", "_")

        # - filter genes and cells
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=10)

        # - add metadata
        adata = add_metadata(adata, name)

        # - save
        adata.write(f"resources/grn_benchmark/inference_data/beeline_{name}.h5ad")
