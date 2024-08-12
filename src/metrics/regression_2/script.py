import pandas as pd
import anndata as ad
import sys
import numpy as np


## VIASH START
par = {
    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
    'layer': 'scgen_pearson',
    'prediction': 'resources/grn-benchmark/negative_control.csv',
    'consensus': 'resources/grn-benchmark/consensus-num-regulators.json',
    'score': 'output/score_regression2.csv',
    'reg_type': 'ridge'
}
## VIASH END

sys.path.append(meta['resources_dir'])
from main import main

if isinstance(par['reg_type'], list) and (len(par['reg_type']) == 1):
    par['reg_type'] = par['reg_type'][0]
assert isinstance(par['reg_type'], str)

print('Reading input data')

output = main(par)

print('Write output to file', flush=True)
print(output)
output.to_csv(par['score'])
print('Completed', flush=True)
