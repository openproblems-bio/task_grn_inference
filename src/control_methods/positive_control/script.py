import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.preprocessing import StandardScaler

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "layer": "lognorm",
  "prior_data":  "resources/grn-benchmark/prior_data.h5ad",
  "output": "resources/grn-benchmark/positive_control.csv",
}
## VIASH END
print('Reading input data')
perturbation_data = ad.read_h5ad(par["perturbation_data"])
gene_names = perturbation_data.var_names.to_numpy()
tf_all = ad.read_h5ad(par["prior_data"]).uns['tf_list']

def create_positive_control(X: np.ndarray) -> np.ndarray:
    X = StandardScaler().fit_transform(X) 
    return np.dot(X.T, X) / X.shape[0]
print('Inferring GRN')
net = create_positive_control(perturbation_data.layers[par["layer"]])

net = pd.DataFrame(net, index=gene_names, columns=gene_names)
net = net.loc[:, net.columns.isin(tf_all)]

pivoted_net = net.reset_index().melt(id_vars='index', var_name='source', value_name='weight')

pivoted_net = pivoted_net.rename(columns={'index': 'target'})
pivoted_net = pivoted_net[pivoted_net['weight'] != 0]
print('Saving')
pivoted_net.to_csv(par["output"])

