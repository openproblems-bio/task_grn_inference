import pandas as pd
import anndata as ad
import sys
import os
import requests
## VIASH START
par = {
  "rna": "resources/grn_benchmark/inference_data/op_rna.h5ad",
  "atac": "resources/grn_benchmark/inference_data/op_atac.h5ad",
  "temp_dir": 'output/scglue_new/',
  "num_workers": 20,
  "prediction": "output/scglue_new/scglue.h5ad",
  "max_n_links": 50000,
  "nes_threshold": 0,
  "rank_threshold": 10000,
  "top_n_targets": 100,
  'normalize': False,
  'extend_range': 150000
}
## VIASH END

import sys



par['annotation_file'] = f"{par['temp_dir']}/gencode.v45.annotation.gtf.gz"
# par['motif_file'] = f"{par['temp_dir']}/JASPAR2022-hg38.bed.gz"
par['motif_file'] = f"{par['temp_dir']}/ENCODE-TF-ChIP-hg38.bed.gz"
# par['motif_file'] = f"output/db/jaspar_encode.bed.gz"

try:
    sys.path.append(meta["resources_dir"])
except:
    meta= {
        "util_dir": 'src/utils/',
        "resources_dir": 'src/methods/multi_omics/scglue'
    }
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["resources_dir"])
from main import main 
from util import parse_args, process_links
par = parse_args(par)


if __name__ == '__main__':
    adata = ad.read_h5ad(par['rna'], backed='r')
    
    print(par)
    net = main(par)

    print('Write output to file', flush=True)
    dataset_id = adata.uns['dataset_id']
    net = net.astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'scglue', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])




