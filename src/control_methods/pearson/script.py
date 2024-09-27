
## VIASH START
par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_0.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'cell_type_specific': False,
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/donor_0_default/pearson.csv',
    "seed": 32,
    'normalize': False
}
## VIASH END
meta = {
    'resources_dir': 'src/utils'
    }

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

import sys
sys.path.append(meta["resources_dir"])
from util import create_corr_net

par['causal'] = True
net = create_corr_net(par)
print('Output GRN')
net.to_csv(par['prediction'])
