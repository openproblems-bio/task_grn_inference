import os
import sys
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

## VIASH START
par = {
    "rna": "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
    "tf_all": "resources_test/prior/tf_all.csv",
    "prediction": "output/geneformer/prediction.h5ad",
    "model": "Geneformer-V2-104M",
    
    "max_n_links": 50000,
    "batch_size": 16,
    "temp_dir": "output/geneformer",
    "file_dir": "src/methods/geneformer",
    "num_genes": 5000,
    "max_cells": 2000,
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
      'util_dir': 'src/utils',
      'helper_dir': 'src/methods/geneformer',
    }
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])
from util import parse_args, efficient_melting
from helper import main

par = parse_args(par)

if __name__ == "__main__":
    net = main(par)
    output = ad.AnnData(
        X=None,
        uns={
            "method_id": "geneformer",
            "dataset_id": dataset_id,
            "prediction": net[["source", "target", "weight", "cell_type"]],
        },
    )
    output.write(par["prediction"])
