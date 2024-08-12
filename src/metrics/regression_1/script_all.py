import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  'max_workers': 4,
  'reg_type': 'ridge',
  'subsample': 2,
}

## VIASH END
print('Reading input data')

meta = {
  "resources_dir":'./'
}
sys.path.append(meta["resources_dir"])
from main import main 
layers = ['pearson', 'lognorm', 'scgen_pearson', 'scgen_lognorm', 'seurat_lognorm', 'seurat_pearson']
grn_models = ['scenicplus','celloracle','figr','granie','scglue','collectri']

os.makedirs('output', exist_ok=True)
for grn_model in grn_models:
  par["score"] = f"output/{grn_model}.csv"
  for ii, layer in enumerate(layers):
    par['prediction'] = f"resources/grn-benchmark/grn_models/{grn_model}.csv"
    par["layer"] = layer

    output = main(par) 
    output.index = [layer]


    if ii == 0:
      score = output
    else:
      score = pd.concat([score,output], axis=0)

    print("Write output to file", flush=True)
    print(grn_model, layer, score)

  print("Write output to file", flush=True)
  score.to_csv(par['score'])