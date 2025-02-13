import pandas as pd
import anndata as ad
import sys
import os
import requests
## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac_d0.h5ad",
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
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--multiomics_atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.multiomics_atac:
    par['multiomics_atac'] = args.multiomics_atac
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
    
# if args.resources_dir:
#     meta['resources_dir'] = args.resources_dir  

# get gene annotation
par['annotation_file'] = f"{par['temp_dir']}/gencode.v45.annotation.gtf.gz"
# par['motif_file'] = f"{par['temp_dir']}/JASPAR2022-hg38.bed.gz"
# par['motif_file'] = f"{par['temp_dir']}/ENCODE-TF-ChIP-hg38.bed.gz"
par['motif_file'] = f"output/db/jaspar_encode.bed.gz"

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
dataset_id = ad.read_h5ad(par['multiomics_rna'], backed='r').uns['dataset_id']
net = net.astype(str)
output = ad.AnnData(X=None, uns={"method_id": 'scglue', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
output.write(par['prediction'])




