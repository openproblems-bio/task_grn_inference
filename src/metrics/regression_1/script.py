import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "prediction": "resources/grn-benchmark/negative_control.csv",
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'lognorm',
  'subsample': 200,
}

## VIASH END
print('Reading input data')

sys.path.append(meta["resources_dir"])
from main import main 

output = main(par) 
print(output)
output.columns = ['S1', 'S2', 'S3', 'S4']

print("Write output to file", flush=True)
print(output)
output.to_csv(par['score'])

# output = ad.AnnData(
#     X=np.empty((0, 0)),
#     uns={
#         "layer": par["layer"],
#         "metric_ids": output.columns.to_numpy(),
#         "metric_values": output.values
#     }
# )
# output.write_h5ad(par["score"], compression="gzip")