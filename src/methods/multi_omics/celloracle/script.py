import pandas as pd
import anndata as ad
import sys
import os

## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
  "temp_dir": 'output/celloracle/',
  "num_workers": 4,
  "prediction": "output/prediction.h5ad",
}
## VIASH END
# meta = {
#   "resources_dir":'resources'
# }
par['links'] = f"{par['temp_dir']}/links.celloracle.links" 

sys.path.append(meta["resources_dir"])
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



