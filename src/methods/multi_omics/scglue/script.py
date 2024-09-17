import pandas as pd
import anndata as ad
import sys
import os

## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna_0.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac_0.h5ad",
  "motif_file": "resources/supplementary/JASPAR2022-hg38.bed.gz",
  "annotation_file": "resources/supplementary/gencode.v45.annotation.gtf.gz",
  "temp_dir": 'output/scglue_ext/',
  "num_workers": 20,
  "prediction": "resources/scglue/prediction.csv",
}
## VIASH END
meta = {
  "resources_dir": 'src/methods/multi_omics/scglue/'
}


sys.path.append(meta["resources_dir"])
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



