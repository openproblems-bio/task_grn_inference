import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
from tqdm import tqdm
import sys
import anndata as ad

## VIASH START
par = {
  'rna': 'resources/grn-benchmark/rna_d0_hvg.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/grnboost2_donor_0_hvg.csv',
  'max_n_links': 50000,
  'cell_type_specific': False,
  'normalize': False,
  'qc': False
}
## VIASH END
if False:
    pass # This block is just a placeholder for the local runs
else:
    sys.path.append(meta["resources_dir"])

from util import process_links, basic_qc

def main(par: dict) -> pd.DataFrame:
    '''
        Main function to infer GRN
    '''
    print('Reading data')
    adata_rna = anndata.read_h5ad(par['rna'])
    if 'qc' in par:
        if par['qc']:
            print('Shape before QC: ', adata_rna.shape)
            adata_rna = basic_qc(adata_rna)
            print('Shape after QC: ', adata_rna.shape)

    gene_names = adata_rna.var_names
    X = adata_rna.X

    # Load list of putative TFs
    tfs = np.loadtxt(par["tf_all"], dtype=str)
    tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

    # GRN inference
    client = Client(processes=False)

    print("Infer grn", flush=True)
  
    network = grnboost2(X, client_or_address=client, gene_names=gene_names, tf_names=tf_names)
    network.rename(columns={'TF': 'source', 'target': 'target', 'importance': 'weight'}, inplace=True)
    network.reset_index(drop=True, inplace=True)
    network = process_links(network, par)    

if __name__ == '__main__':
    net = main(par)   
    # Save inferred GRN
    print('Output GRN')
    # convert the predictions to the benchmark format
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": par['method_id'], "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
