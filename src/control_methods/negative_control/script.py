import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "prediction": "resources/grn_models/default/negative_control.csv",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000
}
## VIASH END
print(par)

def process_links(net, par):
    net = net[net.source!=net.target]
    net = net.sample(par['max_n_links'])
    print(net)
    return net

print('Reading input data')
perturbation_data = ad.read_h5ad(par["perturbation_data"])
gene_names = perturbation_data.var_names.to_numpy()
tf_all = np.loadtxt(par['tf_all'], dtype=str)

n_tf = 500
tfs = tf_all[:n_tf]

def create_negative_control(gene_names) -> np.ndarray:
    ratio = [.98, .01, 0.01]
    
    net = np.random.choice([0, -1, 1], size=((len(gene_names), n_tf)),p=ratio)
    net = pd.DataFrame(net, index=gene_names, columns=tfs)
    return net
print('Inferring GRN')
net = create_negative_control(gene_names)

pivoted_net = net.reset_index().melt(id_vars='index', var_name='source', value_name='weight')

pivoted_net = pivoted_net.rename(columns={'index': 'target'})
pivoted_net = pivoted_net[pivoted_net['weight'] != 0]

pivoted_net = process_links(pivoted_net, par)

print('Saving')
pivoted_net.to_csv(par["prediction"])

