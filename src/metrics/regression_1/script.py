import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "layer": "lognorm",
  "prediction": "resources/grn-benchmark/negative_control.csv",
  'score': r'output/score.csv'
}

## VIASH END
print('Reading input data')

sys.path.append(meta["resources_dir"])
from main import main 

output = main(par) 
output = {key:[value] for key, value in output.items()}

print("Write output to file", flush=True)
output = pd.DataFrame.from_dict(output)
output.to_csv(par['score'])
print("Completed", flush=True)