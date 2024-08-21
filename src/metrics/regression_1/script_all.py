import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler 

par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  'max_workers': 40,
  'reg_type': 'ridge',
  'subsample': -2,
  "tf_all":  "./resources/prior/tf_all.csv",
  "temp_dir": "output"
}

def create_positive_control(X: np.ndarray, groups: np.ndarray):
    grns = []
    for group in tqdm(np.unique(groups), desc="Processing groups"):
        X_sub = X[groups == group, :]
        X_sub = StandardScaler().fit_transform(X_sub)
        grn = np.dot(X_sub.T, X_sub) / X_sub.shape[0]
        grns.append(grn)
    return np.mean(grns, axis=0)
def create_negative_control(gene_names) -> np.ndarray:
    ratio = [.98, .01, 0.01]
    n_tf = 400
    net = np.random.choice([0, -1, 1], size=((len(gene_names), n_tf)),p=ratio)
    net = pd.DataFrame(net, index=gene_names)
    return net


print('Reading input data')
perturbation_data = ad.read_h5ad(par["perturbation_data"])
gene_names = perturbation_data.var_names.to_numpy()
tf_all = np.loadtxt(par['tf_all'], dtype=str)
groups = perturbation_data.obs.cell_type

meta = {
  "resources_dir":'./'
}
sys.path.append(meta["resources_dir"])
from main import main 
layers = ['pearson', 'lognorm', 'scgen_pearson', 'scgen_lognorm', 'seurat_lognorm', 'seurat_pearson']
grn_models = ['scenicplus', 'celloracle', 'figr', 'granie', 'scglue', 'collectri']
controls = ['negative_control', 'positive_control']

os.makedirs(par['temp_dir'], exist_ok=True)
for grn_model in controls + grn_models :
  par["score"] = f"{par['temp_dir']}/reg1-{grn_model}.csv"
  for ii, layer in enumerate(layers):
    par["layer"] = layer
    if grn_model=='positive_control':
      print('Inferring positive control')
      net = create_positive_control(perturbation_data.layers[par["layer"]], groups)

      net = pd.DataFrame(net, index=gene_names, columns=gene_names)
      net = net.loc[:, net.columns.isin(tf_all)]

      pivoted_net = net.reset_index().melt(id_vars='index', var_name='source', value_name='weight')

      pivoted_net = pivoted_net.rename(columns={'index': 'target'})
      pivoted_net = pivoted_net[pivoted_net['weight'] != 0]
      par['prediction'] = f"{par['temp_dir']}/{layer}_positive_control.csv"
      print(par['prediction'])
      pivoted_net.to_csv(par['prediction'])
    elif grn_model=='negative_control':
      print('Inferring negative control')
      net = create_negative_control(gene_names)

      pivoted_net = net.reset_index().melt(id_vars='index', var_name='source', value_name='weight')

      pivoted_net = pivoted_net.rename(columns={'index': 'target'})
      pivoted_net = pivoted_net[pivoted_net['weight'] != 0]
      par['prediction'] = f"{par['temp_dir']}/negative_control.csv"
      pivoted_net.to_csv(par['prediction'])
    else:
      par['prediction'] = f"resources/grn_models/{grn_model}.csv"
    # output = main(par) 
    # output.index = [layer]

    # if ii == 0:
    #   score = output
    # else:
    #   score = pd.concat([score, output], axis=0)

    # print("Write output to file", flush=True)
    # print(grn_model, layer, score)

  # print("Write output to file", flush=True)
  # score.to_csv(par['score'])

