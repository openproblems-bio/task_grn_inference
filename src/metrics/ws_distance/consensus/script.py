import pandas as pd 
import anndata as ad 
from tqdm import tqdm
import numpy as np
import os
import sys


import argparse


arg = argparse.ArgumentParser(description='Compute consensus number of putative regulators for each gene')
arg.add_argument('--dataset', type=str, help='Dataset to use for the analysis')
arg.add_argument('--models_dir', type=str, help='Directory containing the GRN models')
arg.add_argument('--evaluation_data_sc', type=str, help='Path to the evaluation data')
arg.add_argument('--ws_consensus', type=str, help='Path to save the consensus regulators')
arg.add_argument('--tf_all', type=str, help='Path to the file containing all transcription factors')
arg.add_argument('--models', nargs='+', help='List of models to use for the analysis')
args = arg.parse_args()

par = args.__dict__

meta = {
    "resources_dir":'src/metrics/ws_distance/consensus/',
    "utils_dir":'src/utils/'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["utils_dir"])
from helper import main
from util import naming_convention




if __name__ == '__main__':
    par['naming_convention'] = naming_convention
    main(par)
