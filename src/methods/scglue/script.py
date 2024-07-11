import pandas as pd
import anndata as ad

## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
    "output": "output/prediction.csv",
}
## VIASH END

print('Reading input files', flush=True)
multiomics_rna = ad.read_h5ad(par["multiomics_rna"])
    
print('Preprocess data', flush=True)
# ... preprocessing ...

print('Train model', flush=True)
# ... train model ...

print('Generate predictions', flush=True)
# ... generate predictions ...

print('Write output to file', flush=True)
output = pd.DataFrame(
  # ... TODO: fill in data ...
  columns=['source', 'target', 'weight']
)
output.to_csv(par["output"])
