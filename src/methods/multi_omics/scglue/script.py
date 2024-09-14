import pandas as pd
import anndata as ad
import sys
import os

## VIASH START
par = {
  "multiomics_rna": "resources_test/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources_test/grn-benchmark/multiomics_atac.h5ad",
  "motif_file": "resources/supplementary/JASPAR2022-hg38.bed.gz",
  "annotation_file": "resources/supplementary/gencode.v45.annotation.gtf.gz",
  "temp_dir": 'output/scglue/',
  "num_workers": 10,
  "prediction": "output/prediction.csv",
}
## VIASH END
meta = {
  "resources_dir": 'resources'
}


sys.path.append(meta["resources_dir"])
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



