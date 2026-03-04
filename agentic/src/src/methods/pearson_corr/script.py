
import sys
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

## VIASH START
par = {
    'rna': 'resources/grn_benchmark/rna_op.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/prediction.h5ad'
}
## VIASH END

# Load data
rna = ad.read_h5ad(par["rna"])
tf_all = np.loadtxt(par["tf_all"], dtype=str)

# Compute Pearson correlation
genes = rna.var_names
correlations = []
for tf in tf_all:
    if tf in genes:
        for gene in genes:
            if tf != gene:
                corr, _ = pearsonr(rna[:, tf].X.flatten(), rna[:, gene].X.flatten())
                correlations.append((tf, gene, corr))

# Create DataFrame
net = pd.DataFrame(correlations, columns=["source", "target", "weight"])

# Save the inferred network
net['weight'] = net['weight'].astype(str)
output = ad.AnnData(
    X=None,
    uns={
        "method_id": "pearson_corr",
        "dataset_id": "dataset_name",
        "prediction": net[["source", "target", "weight"]]
    }
)
output.write(par["prediction"])
