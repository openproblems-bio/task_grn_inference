
import os 
import anndata as ad 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import scanpy as sc
from sklearn.model_selection import train_test_split 

from scipy.sparse import csr_matrix

## VIASH START
par = {
    'norman_raw': f'resources/datasets_raw/norman.h5ad',
    'norman_bulk': f'resources/extended_data/norman_bulk.h5ad',
    'norman_test_sc': f'resources/grn_benchmark/evaluation_data/norman_sc.h5ad',
    'norman_test_bulk': f'resources/grn_benchmark/evaluation_data/norman_bulk.h5ad',
    'norman_train_sc': f'resources/grn_benchmark/inference_data/norman_rna.h5ad'
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        'resources_dir': 'src/utils/'
    }
    sys.path.append(meta["resources_dir"])

from util import sum_by

def add_metadata(adata):
    adata.uns['dataset_summary'] = 'Single cell RNA-seq data with 231 perturbations (activation) on K562 cells.'
    adata.uns['dataset_description'] = 'How cellular and organismal complexity emerges from combinatorial expression of genes is a central question in biology. High-content phenotyping approaches such as Perturb-seq (single-cell RNA-seq pooled CRISPR screens) present an opportunity for exploring such genetic interactions (GIs) at scale. Here, we present an analytical framework for interpreting high-dimensional landscapes of cell states (manifolds) constructed from transcriptional phenotypes. We applied this approach to Perturb-seq profiling of strong GIs mined from a growth-based, gain-of-function GI map. Exploration of this manifold enabled ordering of regulatory pathways, principled classification of GIs (e.g. identifying suppressors), and mechanistic elucidation of synergistic interactions, including an unexpected synergy between CBL and CNN1 driving erythroid differentiation. Finally, we apply recommender system machine learning to predict interactions, facilitating exploration of vastly larger GI manifolds.'
    adata.uns['data_reference'] = "@article{norman2019exploring, \n\ttitle={Exploring genetic interaction manifolds constructed from rich single-cell phenotypes}, \n\tauthor={Norman, Thomas M and Horlbeck, Max A and Replogle, Joseph M and Ge, Alex Y and Xu, Albert and Jost, Marco and Gilbert, Luke A and Weissman, Jonathan S},\n\tjournal={Science},\n\tvolume={365}, \n\tnumber={6455},\n\tpages={786--793},\n\tyear={2019},\n\tpublisher={American Association for the Advancement of Science}\n\t}"
    adata.uns['data_url'] = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6746554/'
    adata.uns['dataset_id'] = 'norman'
    adata.uns['dataset_name'] = 'Norman'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'X_norm'
    return adata

if __name__ == '__main__':
    # - get the data
    adata = ad.read_h5ad(par['norman_raw'])
    adata = adata[:, ~adata.var_names.duplicated()]
    adata.var_names_make_unique()
    duplicates = adata.var_names[adata.var_names.duplicated()].unique()
    adata = adata[:, ~adata.var_names.isin(duplicates)]
    # - clearn up 
    del adata.obsp 
    del adata.varm
    del adata.uns
    del adata.obsm
    if 'gene_name' in adata.var.columns:
        adata.var = adata.var[['gene_name']]
        adata.var = adata.var.set_index('gene_name')
    else:
        adata.var = adata.var[[]]
    adata.var.index = adata.var.index.astype(str)
    adata.obs = adata.obs[['perturbation', 'is_control', 'perturbation_type']]

    # preprocess  
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=10)

    # - 
    adata.layers['X_norm'] = adata.X.copy()

    # - split to inference and evaluation datasets
    ctr_pertb = adata[adata.obs['is_control']].obs['perturbation'].unique()
    non_ctr_pertubs =adata[~adata.obs['is_control']].obs['perturbation'].unique()
    train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.5, random_state=32)
    train_perturbs = np.concatenate([train_perturbs, ctr_pertb]) # add control perturbations to test set for ws_distance
    test_perturbs = np.concatenate([test_perturbs, ctr_pertb]) 

    adata_train_sc = adata[adata.obs['perturbation'].isin(train_perturbs)] 
    adata_test_sc = adata[adata.obs['perturbation'].isin(test_perturbs)] 


    # - filter genes and cells 
    sc.pp.filter_cells(adata_train_sc, min_genes=100)
    sc.pp.filter_genes(adata_train_sc, min_cells=10)

    sc.pp.filter_cells(adata_test_sc, min_genes=100)
    sc.pp.filter_genes(adata_test_sc, min_cells=10)

    # - pseudo bulk
    adata_bulk = sum_by(adata, unique_mapping=True, col='perturbation') 
    norman_test_bulk = sum_by(adata_test_sc, unique_mapping=True, col='perturbation') # summing over X_norm 

    # - normalize evaluation data
    norman_test_bulk.layers['X_norm'] = norman_test_bulk.X.copy()
    adata_train_sc.layers['X_norm'] = adata_train_sc.X.copy()
    adata_bulk.layers['X_norm'] = adata_bulk.X.copy()

    # - add metadata
    adata_bulk = add_metadata(adata_bulk)
    norman_test_bulk = add_metadata(norman_test_bulk)
    adata_test_sc = add_metadata(adata_test_sc)
    adata_train_sc = add_metadata(adata_train_sc)
    # - save 
    print('saving...')
    adata_bulk.write(par['norman_bulk'])
    adata_test_sc.write(par['norman_test_sc'])
    norman_test_bulk.write(par['norman_test_bulk'])
    adata_train_sc.write(par['norman_train_sc'])