import os
import json

import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse


arg = argparse.ArgumentParser(description='Compute consensus number of putative regulators for each gene')
arg.add_argument('--dataset', type=str, help='Dataset to use for the analysis')
arg.add_argument('--evaluation_data', type=str, help='Path to the evaluation data')
arg.add_argument('--regulators_consensus', type=str, help='Path to save the consensus regulators')
arg.add_argument('--predictions', nargs='+', help='List of models to use for the analysis')
args = arg.parse_args()

par = args.__dict__

meta = {
    "resources_dir":'src/metrics/regression_2/consensus/',
    "utils_dir":'src/utils/'

}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["utils_dir"])

from helper import main

if __name__ == '__main__':
    main(par)

