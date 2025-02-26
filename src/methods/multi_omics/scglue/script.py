import pandas as pd
import anndata as ad
import sys
import os
import requests
## VIASH START
par = {
  "rna": "resources/grn-benchmark/op_rna.h5ad",
  "atac": "resources/grn-benchmark/op_atac.h5ad",
  "temp_dir": 'output/scglue/',
  "num_workers": 20,
  "prediction": "output/scglue_d0_hvg.csv",
  "max_n_links": 50000,
  "nes_threshold": 0,
  "rank_threshold": 10000,
  "top_n_targets": 100,
  'normalize': False,
  'extend_range': 150000
}
## VIASH END

import sys

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
args = parser.parse_args()

for key, value in vars(args).items():
    if value:
        par[key] = value

if __name__ == '__main__':
    adata = ad.read_h5ad(par['rna'], backed='r')
    
    if len(adata) == 2000: # for testing purposes, we will not run the whole pipeline, just keep the format
         net = pd.DataFrame({
            'source': ['A', 'B'],
            'target': ['C', 'D'],
            'weight': [0.1, 0.2],
            'cell_type': ['T cells', 'B cells']})
    else:
        # get gene annotation
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
        print(par)
        net = main(par)

        print('Write output to file', flush=True)
        dataset_id = adata.uns['dataset_id']
    net = net.astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'scglue', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])




