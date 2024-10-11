
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
parser.add_argument('--no_tf_subsetting', action='store_true', default=False, help='Whether to subset based on tf')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--corr_t', type=str, help='Threshold cuttoff for correlation.')
parser.add_argument('--layer', type=str, help='Which layer of adata to use.')


args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.no_tf_subsetting:
    par['no_tf_subsetting'] = args.no_tf_subsetting
if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)
if args.corr_t:
    par['corr_t'] = float(args.corr_t)
if args.layer:
    par['layer'] = args.layer
    
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir  

import sys
sys.path.append(meta["resources_dir"])
from util import create_corr_net

par['causal'] = True
par['normalize'] = True
net = create_corr_net(par)
print('Output GRN')
net.to_csv(par['prediction'])
