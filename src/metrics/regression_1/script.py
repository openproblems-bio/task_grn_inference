import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  # "prediction": "output/pearson_net.csv",
  "prediction": "resources/grn_models/full/portia.csv",
  "method_id": "scenic",
  "min_tf": False,
  "max_n_links": 50000,
  "apply_tf": "true",
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'pearson',
  'subsample': -1,
  'num_workers': 4,
  'skeleton': 'resources/prior/skeleton.csv',
  'apply_skeleton': False,
  'verbose': 0,
  'binarize': True
}
## VIASH END

try:
  sys.path.append(meta["resources_dir"])
except:
  meta = {
    "resources_dir":'src/metrics/regression_1/',
    "util_dir":'src/utils'
  }
  sys.path.append(meta["resources_dir"])
  sys.path.append(meta["util_dir"])

from main import main 

output = main(par) 
print(output)

metric_ids = output.columns.to_numpy()
metric_values = output.values[0]

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

output.write_h5ad(par["score"], compression="gzip")