import pandas as pd
import anndata as ad
import sys

## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
  "motif_file": "resources/grn-benchmark/supp/JASPAR2022-hg38.bed.gz",
  "annotation_file": "resources/grn-benchmark/supp/gencode.v45.annotation.gtf.gz",
  "temp_dir": 'output/scglue/',
  "num_workers": 4,
  "prediction": "output/prediction.csv",
}
## VIASH END
meta = {
  "resources_dir":'resources'
}

sys.path.append(meta["resources_dir"])
from main import main 
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



