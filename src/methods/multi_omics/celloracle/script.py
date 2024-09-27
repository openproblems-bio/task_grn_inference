import pandas as pd
import anndata as ad
import sys
import os
import argparse

## VIASH START
par = {
  "multiomics_rna": "resources_test/grn-benchmark/multiomics_rna_d0_hvg.h5ad",
  "multiomics_atac": "resources_test/grn-benchmark/multiomics_atac.h5ad",
  "base_grn": 'output/celloracle/base_grn.csv',
  "temp_dir": 'output/celloracle/',
  "num_workers": 10,
  "prediction": "output/celloracle_test.h5ad",
}
## VIASH END

parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--multiomics_atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.multiomics_atac:
    par['multiomics_atac'] = args.multiomics_atac
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
    
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir   


par['links'] = f"{par['temp_dir']}/links.celloracle.links" 


import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
    
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir   

sys.path.append(meta["resources_dir"])
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



