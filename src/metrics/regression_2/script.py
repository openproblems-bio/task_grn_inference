import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
  'prediction': 'resources/grn-benchmark/negative_control.csv',
  'score': 'output/score.csv',
  'reg_type': 'ridge'
}

## VIASH END
print('Reading input data')

sys.path.append(meta['resources_dir'])
from main import main 

output = main(par) 
output = output.mean(axis=0).to_frame().T # average across datasets
print(output)
output.columns = [f'S{i + 1}' for i in range(len(output.columns))]
output['Overall'] = output.mean(axis=1)

print('Write output to file', flush=True)
print(output)
output.to_csv(par['score'])
print('Completed', flush=True)
