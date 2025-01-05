import pandas as pd
import anndata as ad
import sys
import numpy as np


## VIASH START
par = {
    'evaluation_data': 'resources/evaluation_datasets/op_perturbation.h5ad',
    'layer': 'X_norm',
    "prediction": "output/models/collectri.csv",
    'tf_all': 'resources/prior/tf_all.csv',
    "max_n_links": 50000,
    'consensus': 'output/models/op_consensus-num-regulators.json',
    'score': 'output/score_regression2.h5ad',
    'reg_type': 'ridge',
    'static_only': True,
    'subsample': -1,
    'num_workers': 4,
    'apply_tf': True,
    'clip_scores': True,
    'method_id': 'grnboost',
    'apply_skeleton': False,
    'skeleton': 'resources/prior/skeleton.csv',
    'verbose': 2
    
}
## VIASH END

print(par)
try:
    sys.path.append(meta['resources_dir'])
except:
    meta = {
    "resources_dir": 'src/metrics/',
    "util": 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util"])

from main import main

print('Reading input data')

output = main(par)

print('Write output to file', flush=True)
print(output)
metric_ids = output.columns.to_numpy()
metric_values = output.values[0]

output = ad.AnnData(
    X=np.empty((0, 0)),
    uns={
        "dataset_id": str(par["dataset_id"]),
        "method_id": f"reg2-{par['method_id']}",
        "metric_ids": metric_ids,
        "metric_values": metric_values
    }
)
output.write_h5ad(par['score'], compression='gzip')
print('Completed', flush=True)
