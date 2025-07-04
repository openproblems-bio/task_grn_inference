import pandas as pd
import anndata as ad
import numpy as np

## VIASH START
par = {
  "rna": "resources/grn-benchmark/rna.h5ad",
  "prediction": "resources/grn_models/default/negative_control.csv",
  "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
  "max_n_links": 50000
}
meta = {
    'resources_dir':'src/control_methods/negative_control'
}
## VIASH END

def process_links(net, par):
    net = net[net.source!=net.target]
    try:
        net = net.sample(par['max_n_links'])
    except:
        net = net.sample(par['max_n_links'], replace=True)
    return net

print('Reading input data')
rna_all = ad.read_h5ad(par["rna_all"])
gene_names = rna_all.var_names.to_numpy()
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

print('Output GRN')
pivoted_net['weight'] = pivoted_net['weight'].astype(str)
output = ad.AnnData(
    X=None,
    uns={
        "method_id": meta['name'],
        "dataset_id": rna_all.uns['dataset_id'],
        "prediction": pivoted_net[["source", "target", "weight"]]
    }
)
output.write_h5ad(par['prediction'])

