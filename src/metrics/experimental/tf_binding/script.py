import os
import sys
import anndata as ad
import numpy as np
import pandas as pd


## VIASH START
par = {
    'prediction': f'resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': f'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'ground_truth': f'resources/grn_benchmark/ground_truth/pbmc.csv',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'genes_n': 5000
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "resources_dir":'src/metrics/experimental/tf_binding/',
    "util_dir": 'src/utils',
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
from helper import main as main 
from util import format_save_score, parse_args

args = parse_args(par)


if __name__ == "__main__":
    output = main(par)
    print(output)

    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    # Custom save for multiple rows
    print('Write output to file', flush=True)
    if not isinstance(output, pd.DataFrame):
        raise ValueError("Expected 'output' to be a pandas DataFrame.")

    # Separate 'gt' column from metric columns
    gt_sources = output['gt'].tolist() if 'gt' in output.columns else []
    metric_cols = [col for col in output.columns if col != 'gt']
    
    metric_ids = metric_cols
    metric_values = output[metric_cols].to_numpy().astype(str)  # Only metric columns

    adata = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": dataset_id,
            "method_id": method_id,
            "metric_ids": metric_ids,
            "metric_values": metric_values,  # store as 2D numpy array of strings
            "gt_sources": gt_sources  # store ground truth sources separately
        }
    )

    adata.write_h5ad(par['score'], compression='gzip')
    print('Completed', flush=True)