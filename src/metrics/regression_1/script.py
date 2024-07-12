import anndata as ad
import pandas as pd
import sys
import os


par = {
  'perturbation_data': r'resources/grn-benchmark/perturbation_data.h5ad',
  'layer': r'lognorm',
  'prediction': r'resources/grn-benchmark/collectri.csv',
  'score': r'output/score.csv',
  'reg_type': r'ridge',
  'subsample': 200
}
meta = {
    "resources_dir": "src/metrics/regression_1",
}
sys.path.append(meta["resources_dir"])

from main import main 

output = main(par) 
output = {key:[value] for key, value in output.items()}

print("Write output to file", flush=True)
output = pd.DataFrame.from_dict(output)
output.to_csv(par['score'])
print("Completed", flush=True)
