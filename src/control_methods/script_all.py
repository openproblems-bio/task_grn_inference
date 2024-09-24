import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 
import random

par = {
  'write_dir': "resources/grn_models/d0_hvgs",
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad",
  
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000,
  'layer': 'scgen_pearson',
  'cell_type_specific': False,
  'normalize': False,
  'causal': True
}

meta = {
  "resources_dir": 'src/control_methods',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['write_dir'], exist_ok=True)

from util import create_corr_net


#---- run for pearson_corr
par['prediction'] = f"{par['write_dir']}/pearson_corr.csv"
par['causal'] = True
net = create_corr_net(par)
net.to_csv(par['prediction'])
#---- run for negative control
from negative_control.main import main 
par['prediction'] = f"{par['write_dir']}/negative_control.csv"
net = main(par)
net.to_csv(par['prediction'])
#---- run for positive_control
par['multiomics_rna'] = par['perturbation_data']
par['prediction'] = f"{par['write_dir']}/positive_control.csv"
net = create_corr_net(par)
net.to_csv(par['prediction'])
