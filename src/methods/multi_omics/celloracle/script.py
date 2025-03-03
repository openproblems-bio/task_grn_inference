import pandas as pd
import anndata as ad
import sys
import os
import argparse

## VIASH START
par = {
  "rna": "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
  "atac": "resources_test/grn_benchmark/inference_data/op_atac.h5ad",
  "base_grn": 'output/celloracle/base_grn.csv',
  "temp_dir": 'output/celloracle/',
  "num_workers": 10,
  "prediction": "output/celloracle.h5ad"}
## VIASH END

parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=int, help='Number of cores')
parser.add_argument('--max_n_links', type=int)
args = parser.parse_args()

for key in vars(args):
    if getattr(args, key) is not None:
        par[key] = getattr(args, key)
    

try:
    sys.path.append(meta["resources_dir"])
except: # for local run
    meta = {'resources_dir': './'}
    sys.path.append(meta["resources_dir"])

from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)

 
if 'base_grn' not in par:
    par['base_grn'] = f"{par['temp_dir']}/base_grn.csv" 
if 'links' not in par:
    par['links'] = f"{par['temp_dir']}/links.celloracle.links"

if __name__ == '__main__':
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
    net = main(par)

    print('Write output to file', flush=True)
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'celloracle', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])



