import os
import anndata
import numpy as np
import pandas as pd
import subprocess
import ast
import requests
import scipy.sparse as sp
import sys
import anndata as ad
import argparse

## VIASH START
par = {
  'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
  "tf_all": 'resources/grn_benchmark/prior/tf_all.csv',
  'prediction': 'output/grnboost2_test.h5ad',
  'max_n_links': 50000,
  'cell_type_specific': False,
  'normalize': False,
  'num_workers': 10,
  'temp_dir': 'output/temdir/',
  'seed': 42,
  'qc': False,
  'layer': 'X_norm'
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--normalize', action='store_true')
args = parser.parse_args()

par_local = vars(args)

for key, value in par_local.items():
    if value is not None:
        par[key] = value

## LOCAL END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    'resources_dir': 'src/utils'
    }
    sys.path.append(meta["resources_dir"])

from util import process_links

def run_grn(par):
    print('Run grn')
    command = [
        "pyscenic", "grn",
        "--num_workers", str(par['num_workers']),
        "--seed", str(par['seed']),
        "-o", str(par['expr_mat_adjacencies']),
        "--method", "grnboost2", 
        str(par['expression_data']),
        par['tf_all']
    ]
    subprocess.run(command, check=True)
def format_data(par):
    print('Read data')
    adata_rna = ad.read_h5ad(par['rna'])  
    gene_names = adata_rna.var_names
    if sp.issparse(adata_rna.X):
        adata_rna.X = adata_rna.X.toarray()
    pd.DataFrame(adata_rna.X, columns=gene_names).to_csv(par['expression_data'], sep='\t', index=False)
  

def main(par):
    os.makedirs(par['temp_dir'], exist_ok=True)
    par['expr_mat_adjacencies'] =  os.path.join(par['temp_dir'], "expr_mat_adjacencies.tsv")
    par['expression_data'] = os.path.join(par['temp_dir'], "expression_data.tsv")

    format_data(par)
    run_grn(par)
    network = pd.read_csv(par['expr_mat_adjacencies'], sep='\t') 
    network.rename(columns={'TF': 'source', 'importance': 'weight'}, inplace=True)
    network.reset_index(drop=True, inplace=True)
    network = process_links(network, par)
    return network


if __name__ == '__main__':
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']

    net = main(par)   

    # Save inferred GRN
    print('Output GRN')
    # convert the predictions to the benchmark format
    net = net.astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'grnboost2', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
