import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 
import random

par = {
  'models_dir': "resources/grn_models/d0_hvgs",
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad",
  
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000,
  'layer': 'scgen_pearson',
  'cell_type_specific': False,
  'normalize': False,
  'causal': True
}

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--models_dir', type=str, help='where to write the model')
args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.models_dir:
    par['models_dir'] = args.models_dir
    
if args.resources_dir:
    meta = {'resources_dir': args.resources_dir} 

meta = {
  "resources_dir": 'src/control_methods',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['models_dir'], exist_ok=True)

from util import create_corr_net

print(par)
#---- run for pearson_corr
print("run for pearson_corr")
par['prediction'] = f"{par['models_dir']}/pearson_corr.csv"
par['causal'] = True
net = create_corr_net(par)
net.to_csv(par['prediction'])
#---- run for negative control
print("run for negative control")
from negative_control.main import main 
par['prediction'] = f"{par['models_dir']}/negative_control.csv"
net = main(par)
net.to_csv(par['prediction'])
#---- run for positive_control
print("run for positive control")
par['multiomics_rna'] = par['perturbation_data']
par['prediction'] = f"{par['models_dir']}/positive_control.csv"
net = create_corr_net(par)
net.to_csv(par['prediction'])
