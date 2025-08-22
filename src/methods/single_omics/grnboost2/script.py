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
  'rna': 'resources/grn_benchmark/inference_data/replogle_rna.h5ad',
  "tf_all": 'resources/grn_benchmark/prior/tf_all.csv',
  'prediction': 'output/grnboost2_test.h5ad',
  'temp_dir': 'output/grnboost2',
  'num_workers': 1,
  'seed': "32",
  'normalize': False
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
      'util_dir': 'src/utils',
      'helper_dir': 'src/methods/single_omics/grnboost2',
    }
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])
from util import process_links
from helper import format_data, run_grn


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
  
if __name__=='__main__':
  dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
  net = main(par)
  
  net['weight'] = net['weight'].astype(str)
  output = ad.AnnData(X=None, uns={"method_id": "grnboost", "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
  output.write(par['prediction'])



