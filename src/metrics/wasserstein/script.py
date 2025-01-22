import pandas as pd
import anndata as ad
import sys
import numpy as np


## VIASH START
par = {
    'prediction': 'resources/grn_models/adamson/pearson_corr.csv',
    'evaluation_data_sc': 'resources/datasets_raw/adamson_sc_counts.h5ad',
    'score': 'output/score.h5ad',
    'ws_consensus': 'resources/prior/consensus_ws_distance_adamson.csv',
    'ws_distance_background':'resources/prior/ws_distance_background_adamson.csv',
    'layer': 'X_norm',
    'dataset_id': 'dataset_id',
    'method_id': 'pearson_corr'
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
    from main import main 
except:
    meta = {
    "resources_dir": 'src/metrics/wasserstein',
    "util": 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util"])
    from main import main 

if __name__ == '__main__':
    _, mean_scores = main(par)
    
    output = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": str(par["dataset_id"]),
            "method_id": f"reg2-{par['method_id']}",
            "metric_ids": mean_scores.columns.values,
            "metric_values": mean_scores.values[0]
        }
    )
    output.write_h5ad(par['score'], compression='gzip')
    print('Completed', flush=True)