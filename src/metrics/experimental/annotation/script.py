import os
import sys
import anndata as ad
import numpy as np
import pandas as pd
import argparse


## VIASH START
par = {
    'prediction': f'resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': f'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'pathway_file': '/vol/projects/jnourisa/prior/h.all.v2024.1.Hs.symbols.gmt',
    'fdr_threshold': 0.05,
    'min_pathway_size': 5,
    'max_pathway_size': 500,
    'min_targets': 10,
    'max_targets': 100,  # Top K edges by absolute weight
    'run_gene_set_recovery': True,  # Enable ULM-based gene set recovery
    'ulm_activity_threshold': 0.0,
    'ulm_pvalue_threshold': 0.01,
    'ulm_baseline_method': 'zero_centered'
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "resources_dir":'src/metrics/experimental/annotation/',
    "util_dir": 'src/utils',
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
from helper import main as main_helper
from util import format_save_score, parse_args

# Custom argument parser for pathway-specific parameters
def parse_pathway_args(par):
    """Parse pathway-specific arguments and add them to par dict"""
    parser = argparse.ArgumentParser()
    
    # Standard arguments (from util.parse_args)
    parser.add_argument('--rna', type=str)
    parser.add_argument('--rna_all', type=str)
    parser.add_argument('--atac', type=str)
    parser.add_argument('--prediction', type=str)
    parser.add_argument('--score', type=str)
    parser.add_argument('--layer', type=str)
    parser.add_argument('--temp_dir', type=str)
    parser.add_argument('--tf_all', type=str)
    parser.add_argument('--skeleton', type=str)
    parser.add_argument('--apply_skeleton', action='store_true')
    parser.add_argument('--apply_tf', action='store_true')
    parser.add_argument('--max_n_links', type=int)
    parser.add_argument('--reg_type', type=str)
    parser.add_argument('--num_workers', type=int)
    parser.add_argument('--regulators_consensus', type=str)
    parser.add_argument('--evaluation_data', type=str)
    parser.add_argument('--ws_consensus', type=str)
    parser.add_argument('--ws_distance_background', type=str)
    parser.add_argument('--group_specific', type=str)
    parser.add_argument('--evaluation_data_de', type=str)
    parser.add_argument('--evaluation_data_sc', type=str)
    parser.add_argument('--ground_truth_unibind', type=str)
    parser.add_argument('--ground_truth_chipatlas', type=str)
    parser.add_argument('--ground_truth_remap', type=str)
    
    # Pathway-specific arguments
    parser.add_argument('--pathway_file', type=str, 
                       default='/vol/projects/jnourisa/prior/h.all.v2024.1.Hs.symbols.gmt')
    parser.add_argument('--fdr_threshold', type=float, default=0.05)
    parser.add_argument('--min_pathway_size', type=int, default=5)
    parser.add_argument('--max_pathway_size', type=int, default=500)
    parser.add_argument('--min_targets', type=int, default=10)
    parser.add_argument('--max_targets', type=int, default=100,
                       help='Maximum targets per TF (top K by absolute edge weight)')
    
    # Gene set recovery arguments (ULM-based)
    parser.add_argument('--run_gene_set_recovery', action='store_true',
                       help='Run ULM-based gene set recovery metric')
    parser.add_argument('--ulm_activity_threshold', type=float, default=0.0,
                       help='Minimum activity score for pathway to be considered active')
    parser.add_argument('--ulm_pvalue_threshold', type=float, default=0.01,
                       help='P-value threshold for pathway activity significance')
    parser.add_argument('--ulm_baseline_method', type=str, default='zero_centered',
                       choices=['zero_centered', 'permutation', 'random_genesets'],
                       help='Method for determining pathway activity baseline')
    
    args = parser.parse_args()
    for k, v in vars(args).items():
        if v is not None:
            par[k] = v
    return par

args = parse_pathway_args(par)


if __name__ == "__main__":
    output = main_helper(par)
    print(output)

    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    format_save_score(output, method_id, dataset_id, par['score'])