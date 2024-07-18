import pandas as pd
import anndata as ad
import sys
import subprocess
import os

## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
  "temp_dir": 'output/temp_figr/',
  "num_workers": 4,
  "n_topics":48,
  "prediction": "output/prediction.csv",
}
## VIASH END
# meta = {
#   "resources_dir":'resources'
# }
print('\n',par)

sys.path.append(meta["resources_dir"])
print("Formating input data")
from format_data import format_data 
# format_data(par)

# Convert the parameters to command line arguments
args = ["Rscript", "format_data_r.R"]
for key, value in par.items():
    args.append(f"--{key}")
    args.append(str(value))

# Run the R script with the command line arguments
result = subprocess.run(args)

# Check if the subprocess was successful
if result.returncode != 0:
    raise RuntimeError(f"Subprocess failed with return code {result.returncode}")

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])
