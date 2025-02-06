import pandas as pd
import anndata as ad
import sys
import os
import argparse

## VIASH START
par = {
  "rna": "resources/grn_benchmark/inference_data/op_rna.h5ad",
  "atac": "resources/grn_benchmark/inference_data/op_atac.h5ad",
  "base_grn": 'output/celloracle/base_grn.csv',
  "temp_dir": 'output/celloracle/',
  "num_workers": 10,
  "prediction": "output/celloracle.h5ad"}
## VIASH END

parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=int, help='Number of cores')
parser.add_argument('--max_n_links', type=int)
args = parser.parse_args()

if args.rna:
    par['rna'] = args.rna
if args.atac:
    par['atac'] = args.atac
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.max_n_links:
    par['max_n_links'] = args.max_n_links
    
if args.resources_dir:
    meta = {}
    meta['resources_dir'] = args.resources_dir   
try:
    sys.path.append(meta["resources_dir"])
except:
    pass
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)

 
if 'base_grn' not in par:
    par['base_grn'] = f"{par['temp_dir']}/base_grn.csv" 
if 'links' not in par:
    par['links'] = f"{par['temp_dir']}/links.celloracle.links"

net = main(par)

print('Write output to file', flush=True)
net['weight'] = net['weight'].astype(str)
output = ad.AnnData(X=None, uns={"method_id": par['method_id'], "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
output.write(par['prediction'])



