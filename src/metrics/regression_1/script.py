import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  "prediction": "resources/grn_models/donor_0_default/pearson_causal.csv",
  "method_id": "scenic",
  "min_tf": False,
  "max_n_links": 50000,
  "apply_tf": "true",
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'scgen_pearson',
  'subsample': -2,
  'num_workers': 4,
}
## VIASH END
# meta = {
#   "resources_dir":'src/metrics/regression_1/'
# }
sys.path.append(meta["resources_dir"])
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