import pandas as pd
import anndata as ad
import sys
import numpy as np
import sys 

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "prediction": "resources/grn_models/default/negative_control.csv",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000
}
## VIASH END
print(par)
meta = {
  'resources_dir':'src/control_methods/negative_control'
}
sys.path.append(meta['resources_dir'])
from main import main
prediction = main(par)
prediction.to_csv(par["prediction"])

