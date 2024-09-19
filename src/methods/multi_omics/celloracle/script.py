import pandas as pd
import anndata as ad
import sys
import os

## VIASH START
par = {
  "multiomics_rna": "resources_test/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources_test/grn-benchmark/multiomics_atac.h5ad",
  "base_grn": 'output/celloracle/base_grn.csv',
  "temp_dir": 'output/celloracle/',
  "num_workers": 4,
  "prediction": "output/celloracle_test.h5ad",
}
## VIASH END
meta = {
  "resources_dir":'src/methods/multi_omics/celloracle'
}
par['links'] = f"{par['temp_dir']}/links.celloracle.links" 

sys.path.append(meta["resources_dir"])
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



