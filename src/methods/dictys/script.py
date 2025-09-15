
import pandas as pd
import numpy as np
import os
import argparse
import warnings
import sys
import importlib

warnings.filterwarnings("ignore")

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
    'temp_dir': 'output/temp/',
    'prediction': 'output/temp/prediction.h5ad',
    'num_workers': 1
}
## VIASH END
try: 
    sys.path.append(meta["resources_dir"])
    par['frag_to_bam'] = f"{meta['resources_dir']}/frag_to_bam.py"

except:
    meta = {
        'utils_dir': 'src/utils',
        'helper_dir': 'src/methods/dictys/',
        'frag_to_bam': 'src/methods/dictys/frag_to_bam.py'
    }
    sys.path.append(meta['utils_dir'])
    sys.path.append(meta['helper_dir'])
    par['frag_to_bam'] = meta['frag_to_bam']

from helper import main 
# from util import parse_args


# par = parse_args(par)

if __name__ == '__main__':
    main(par)