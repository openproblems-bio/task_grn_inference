import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "prediction": "resources/grn-benchmark/grn_models/collectri.csv",
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'lognorm',
  'subsample': 200,
  'max_workers': 4,
}

## VIASH END
sys.path.append(meta["resources_dir"])
from main import main 

output = main(par) 
print(output)
# output.columns = ['S1', 'S2', 'S3', 'S4']
# output.index=[par["layer"]]
# print("Write output to file", flush=True)
# print(output)
# output.to_csv(par['score'])


metric_ids = output.columns.to_numpy()
metric_values = output.values[0]
# if metric_ids.ndim == 1:
#     metric_ids = metric_ids.reshape(1, -1)
# if metric_values.ndim == 1:
#     metric_values = metric_values.reshape(1, -1)

print(metric_ids.shape, metric_values.shape)
output = ad.AnnData(
    X=np.empty((0, 0)),
    uns={
        "dataset_id": par["layer"],
        "method_id": f"reg1-{par['method_id']}",
        "metric_ids": metric_ids,
        "metric_values": metric_values
    }
)


# print(output.uns)

output.write_h5ad(par["score"], compression="gzip")