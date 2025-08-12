import pandas as pd
import anndata as ad
import sys
import numpy as np
import random 
import os
import warnings
warnings.simplefilter("ignore")

## VIASH START
par = {

}
## VIASH END

print(par)

if __name__ == '__main__':
    
    results = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": 'dummy_dataset',
            "method_id": 'dummy_method',
            "metric_ids": ['dummy_metric'],
            "metric_values": [0]
        }
    )
    print(results.uns)

    results.write_h5ad(par["score"], compression="gzip")