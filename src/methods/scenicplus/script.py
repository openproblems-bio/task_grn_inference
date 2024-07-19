import pandas as pd
import anndata as ad
import sys

## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
  "prediction": "output/prediction.csv",
}
## VIASH END
sys.path.append(meta["resources_dir"])
from main import main 
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



