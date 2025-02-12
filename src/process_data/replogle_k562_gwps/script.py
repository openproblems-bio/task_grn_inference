import sys
import anndata as ad

import argparse
import numpy as np
from sklearn.model_selection import train_test_split 
import scanpy as sc
import gc

## VIASH START
argparser = argparse.ArgumentParser()
argparser.add_argument('--input', type=str, required=True)
argparser.add_argument('--tf_all', type=str, required=True)
argparser.add_argument('--adata_test_sc', type=str, required=True)
argparser.add_argument('--adata_train_sc', type=str, required=True)
argparser.add_argument('--adata_train_sc_subset', type=str, required=True)

args = argparser.parse_args()
par = vars(args)

par_local = {
    'test_perturbs': f'resources/grn_benchmark/prior/replogle_test_perturbs.csv',

}
## VIASH END

meta = {
    'resources_dir': 'src/utils/'
}
sys.path.append(meta["resources_dir"])

from util import sum_by, fetch_gene_info


def format_raw_data(adata: ad.AnnData) -> ad.AnnData:
    '''
    Format the raw data
    '''
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    adata.var['gene_name'] = adata.var['gene_name'].astype(str)
    adata.var = adata.var.reset_index().set_index('gene_name')[['gene_id','chr', 'start', 'end']]
    adata.obs = adata.obs.rename(columns={'gene': 'perturbation'})[['perturbation']]
    adata.obs['perturbation_type'] = 'knockdown'
    adata.obs['is_tf'] = adata.obs['perturbation'].isin(tf_all)
    adata.var_names_make_unique()  
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    adata.obs['is_control'] = adata.obs['perturbation']=='non-targeting'
    
    return adata

def split_data(adata: ad.AnnData):
    '''
    Split the data into train and test
    
    '''
    unique_perts = adata.obs['perturbation'].unique()
    
    tf_all = adata.obs[adata.obs['is_tf']]['perturbation'].unique()
    test_perturbs = np.loadtxt(par_local['test_perturbs'], dtype=str)
    test_tfs = np.intersect1d(tf_all, test_perturbs)
    train_pertubs = np.setdiff1d(unique_perts, test_perturbs)

    ctrs =  ['non-targeting']
    test_tfs = np.concatenate([test_tfs, ctrs])
    
    print(f"Test TFs: {len(test_tfs)}")

    adata_test = adata[adata.obs['perturbation'].isin(test_tfs)]
    adata_train = adata[adata.obs['perturbation'].isin(train_pertubs)] # all perturbs

    adata_train_short = adata[adata.obs['perturbation'].isin(train_pertubs[:500])]

    print(f"Train : {adata_train.shape}, Test: {adata_test.shape}, Train_short: {adata_train_short.shape}")

    return adata_train, adata_test, adata_train_short



def normalize(adata: ad.AnnData) -> ad.AnnData:
    X_norm = sc.pp.normalize_total(adata, inplace=False)['X']
    X_norm = sc.pp.log1p(X_norm, copy=True)

    adata.layers['X_norm'] = X_norm
    return adata
def normalize_pearson(adata: ad.AnnData) -> ad.AnnData:
    X_norm = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)['X']

    adata.layers['X_norm'] = X_norm
    return adata


def main(par):
    adata = ad.read_h5ad(par['input'], backed='r') 
    if False: # to test 
        tf_all = np.loadtxt(par['tf_all'], dtype=str)
        filter = adata.obs['gene'][adata.obs['gene'].isin(tf_all)].unique()[0:10]
        adata = adata[adata.obs['gene'].isin(filter)].to_memory()  

    # Format the raw data
    print('Formatting raw data...')
    adata = format_raw_data(adata)

    # - process adata
    print('Processing data...')
    adata_sc = adata.to_memory()
    if False:
        adata_bulk = psedudobulk_fun(adata_sc)
        adata_bulk = normalize(adata_bulk)
        adata_bulk.write(par['adata_bulk'])
        del adata_bulk, adata_sc 
        gc.collect()
    # - data split
    print('Splitting data...')
    adata_train_sc, adata_test_sc, adata_train_sc_subset = split_data(adata)
    del adata
    gc.collect()

    # - process inference data
    print('Process inference data...')
    adata_train_sc = adata_train_sc.to_memory()
    adata_train_sc = normalize(adata_train_sc)
    if False:
        adata_train_bulk = sum_by(adata_train_sc, col='perturbation', unique_mapping=True)
        adata_train_bulk = normalize_pearson(adata_train_bulk)
        adata_train_bulk.uns['dataset_id'] = 'replogle'
        adata_train_bulk.write(par['adata_train_bulk'])
        
        
        del adata_train_bulk
        gc.collect()

    if False:
        print('adata_train_bulk: ', adata_train_bulk.shape)
        print('pertrbations in adata_train_bulk: ', adata_train_bulk.obs['perturbation'].nunique())
        print('adata_train_sc: ', adata_train_sc.shape)
        print('pertrbations in adata_train_sc: ', adata_train_sc.obs['perturbation'].nunique())
    
    adata_train_sc.uns['dataset_id'] = 'replogle'
    adata_train_sc.write(par['adata_train_sc'])

    adata_train_sc_subset.uns['dataset_id'] = 'replogle'
    adata_train_sc_subset.write(par['adata_train_sc_subset'])
    del adata_train_sc, adata_train_sc_subset
    gc.collect()

    # - process test data
    print('Process test data...')
    adata_test_sc = adata_test_sc.to_memory()
    adata_test_sc = normalize(adata_test_sc)
    adata_test_sc.uns['dataset_id'] = 'replogle'
    adata_test_sc.write(par['adata_test_sc'])

    if False:
        adata_test_bulk = sum_by(adata_test_sc, col='perturbation', unique_mapping=True)
        adata_test_bulk.uns['dataset_id'] = 'replogle'
        adata_test_bulk = normalize_pearson(adata_test_bulk)
        adata_test_bulk.write(par['adata_test_bulk'])


if __name__ == '__main__':
    main(par)

