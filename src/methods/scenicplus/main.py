import os 
import pandas as pd 
import anndata as ad
import scenicplus

def main(par):
    temp_dir = 'output/temp_scenicplus/'
    num_workers = 4
    os.makedirs(temp_dir, exist_ok=True)
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['multiomics_rna'])
    atac = ad.read_h5ad(par['multiomics_atac'])
    print('Preprocess data', flush=True)
    
    grn = pd.DataFrame(
    data = {'source':['tf1'], 'target':['g1'], 'weight':[1]}
    )

    return grn
    