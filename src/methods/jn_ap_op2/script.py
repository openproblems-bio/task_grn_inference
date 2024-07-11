# This script is based an IPython notebook:
# https://github.com/AntoinePassemiers/Open-Challenges-Single-Cell-Perturbations/blob/master/op2-de-dl.ipynb

import sys
import pandas as pd
import anndata as ad
import os
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'  # Make PyTorch deterministic on GPU
import numpy as np
import pandas as pd
import tqdm
import torch

## VIASH START
par = {
    "de_train": "resources/neurips-2023-kaggle/de_train.h5ad",
    "id_map": "resources/neurips-2023-kaggle/id_map.csv",
    "output": "output.h5ad",
    "n_replica": 1,
    "submission_names": ["dl40"]
}
meta = {
    "resources_dir": "src/methods/jn_ap_op2",
}
## VIASH END

sys.path.append(meta["resources_dir"])

from helper import plant_seed, MultiOutputTargetEncoder, train

print('Reading input files', flush=True)
de_train_h5ad = ad.read_h5ad(par["de_train_h5ad"])
id_map = pd.read_csv(par["id_map"])

gene_names = list(de_train_h5ad.var_names)

print('Preprocess data', flush=True)
SEED = 0xCAFE
USE_GPU = True
if USE_GPU and torch.cuda.is_available():
    print('using device: cuda')
else:
    print('using device: cpu')
    USE_GPU = False
    
# Make Python deterministic?
os.environ['PYTHONHASHSEED'] = str(int(SEED))

# Make PyTorch deterministic
torch.use_deterministic_algorithms(True)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.enabled = False
torch.set_num_threads(1)
plant_seed(SEED, USE_GPU)

print('Data location', flush=True)
# Data location
cell_types = de_train_h5ad.obs['cell_type'].astype(str)
sm_names = de_train_h5ad.obs['sm_name'].astype(str)

data = de_train_h5ad.layers[par["layer"]]

print('Train model', flush=True)
# ... train model ...
encoder = MultiOutputTargetEncoder()

encoder.fit(np.asarray([cell_types, sm_names]).T, data)

X = torch.FloatTensor(encoder.transform(np.asarray([cell_types, sm_names]).T))
X_submit = torch.FloatTensor(encoder.transform(np.asarray([id_map.cell_type, id_map.sm_name]).T))

if USE_GPU:
    X = X.cuda()

print('Generate predictions', flush=True)
# ... generate predictions ...

Y_submit_ensemble = []
for SUBMISSION_NAME in par["submission_names"]:
  #train the models and store them
  models = []
  for i in range(par["n_replica"] ):
      seed = i
      if SUBMISSION_NAME == 'dl40':
          model = train(X, torch.FloatTensor(data), np.arange(len(X)), seed, n_iter=40, USE_GPU=USE_GPU)
      elif SUBMISSION_NAME == 'dl200':
          model = train(X, torch.FloatTensor(data), np.arange(len(X)), seed, n_iter=200, USE_GPU=USE_GPU)
      else:
          model = train(X, torch.FloatTensor(data), np.arange(len(X)), seed, n_iter=40, USE_GPU=USE_GPU)
      model.eval()
      models.append(model)
      torch.cuda.empty_cache()
  # predict 
  Y_submit =  []
  for i, x in tqdm.tqdm(enumerate(X_submit), desc='Submission'):
    # Predict on test sample using a simple ensembling strategy:
    # take the median of the predictions across the different models
    y_hat = []
    for model in models:
        model = model.cpu()
        y_hat.append(np.squeeze(model.forward(x.unsqueeze(0)).cpu().data.numpy()))
    y_hat = np.median(y_hat, axis=0)

    values = [f'{x:.5f}' for x in y_hat]
    Y_submit.append(values)
    
  Y_submit_ensemble.append(np.asarray(Y_submit).astype(np.float32))
    
Y_submit_final = np.mean(Y_submit_ensemble, axis=0)

print('Write output to file', flush=True)
output = ad.AnnData(
    layers={"prediction": Y_submit_final},
    obs=pd.DataFrame(index=id_map["id"]),
    var=pd.DataFrame(index=gene_names),
    uns={
      "dataset_id": de_train_h5ad.uns["dataset_id"],
      "method_id": meta["functionality_name"]
    }
)

output.write_h5ad(par["output"], compression="gzip")
