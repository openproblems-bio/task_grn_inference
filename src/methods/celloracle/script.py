import pandas as pd
import anndata as ad
import sys
import os

## VIASH START
par = {
  "rna": "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
  "atac": "resources_test/grn_benchmark/inference_data/op_atac.h5ad",
  "base_grn": 'output/celloracle/base_grn.csv',
  "annotated_peaks": 'output/celloracle/annotated_peaks.csv',
  "temp_dir": 'output/celloracle/',
  "num_workers": 10,
  "prediction": "output/celloracle.h5ad"}
## VIASH END
import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--rna', type=str, help='Path to the input RNA data in h5ad format.')
argparser.add_argument('--atac', type=str, help='Path to the input ATAC data in h5ad format.')
argparser.add_argument('--prediction', type=str, help='Path to the output prediction in h5ad format.')
argparser.add_argument('--annotated_peaks', type=str, default=par['annotated_peaks'], help='Path to store the annotated peaks.')
args = argparser.parse_args()
if args.rna is not None:
    par['rna'] = args.rna
if args.atac is not None:
    par['atac'] = args.atac
if args.prediction is not None:
    par['prediction'] = args.prediction


try:
    sys.path.append(meta["resources_dir"])
except: # for local run
    meta = {'resources_dir': './'}
    sys.path.append(meta["resources_dir"])

from helper import main 
os.makedirs(par['temp_dir'], exist_ok=True)

 
if 'base_grn' not in par:
    par['base_grn'] = f"{par['temp_dir']}/base_grn.csv" 
if 'links' not in par:
    par['links'] = f"{par['temp_dir']}/links.celloracle.links"

if __name__ == '__main__':
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
    net = main(par)

    print('Write output to file', flush=True)
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'celloracle', "dataset_id": dataset_id, "prediction": net})
    output.write(par['prediction'])



