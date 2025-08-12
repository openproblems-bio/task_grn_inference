
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
    'adamson_raw': f'resources/datasets_raw/adamson.h5ad',
    'adamson_bulk': f'resources/extended_data/adamson_bulk.h5ad',
    'adamson_test_sc': f'resources/grn_benchmark/evaluation_data/adamson_sc.h5ad',
    'adamson_test_bulk': f'resources/grn_benchmark/evaluation_data/adamson_bulk.h5ad',
    'adamson_train_sc': f'resources/grn_benchmark/inference_data/adamson_rna.h5ad',
    'adamson_train_sc_test': f'resources_test/grn_benchmark/inference_data/adamson_rna.h5ad',
    'adamson_test_bulk_test': f'resources_test/grn_benchmark/evaluation_data/adamson_bulk.h5ad',
    'adamson_test_sc_test': f'resources_test/grn_benchmark/evaluation_data/adamson_sc.h5ad'

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
    adata.uns['dataset_summary'] = 'Single cell RNA-seq data with 82 perturbations (KD) on K562 cells.'
    adata.uns['dataset_description'] = 'Functional genomics efforts face tradeoffs between number of perturbations examined and complexity of phenotypes measured. We bridge this gap with Perturb-seq, which combines droplet-based single-cell RNA-seq with a strategy for barcoding CRISPR-mediated perturbations, allowing many perturbations to be profiled in pooled format. We applied Perturb-seq to dissect the mammalian unfolded protein response (UPR) using single and combinatorial CRISPR perturbations. Two genome-scale CRISPR interference (CRISPRi) screens identified genes whose repression perturbs ER homeostasis. Subjecting ∼100 hits to Perturb-seq enabled high-precision functional clustering of genes. Single-cell analyses decoupled the three UPR branches, revealed bifurcated UPR branch activation among cells subject to the same perturbation, and uncovered differential activation of the branches across hits, including an isolated feedback loop between the translocon and IRE1α. These studies provide insight into how the three sensors of ER homeostasis monitor distinct types of stress and highlight the ability of Perturb-seq to dissect complex cellular responses.'
    adata.uns['data_reference'] = "@article{adamson2016multiplexed,\n\ttitle={A multiplexed single-cell CRISPR screening platform enables systematic dissection of the unfolded protein response},\n\tauthor={Adamson, Britt and Norman, Thomas M and Jost, Marco and Cho, Min Y and Nuez, James K and Chen, Yuwen and Villalta, Jacqueline E and Gilbert, Luke A and Horlbeck, Max A and Hein, Marco Y and others},\n\tjournal={Cell},\n\tvolume={167},\n\tnumber={7},\n\tpages={1867--1882},\n\tyear={2016},\n\tpublisher={Elsevier}}"
    adata.uns['data_url'] = 'https://pubmed.ncbi.nlm.nih.gov/27984733/'
    adata.uns['dataset_id'] = 'adamson'
    adata.uns['dataset_name'] = 'Adamson'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    return adata

def main(par):
    # - get the data
    adata = ad.read_h5ad(par['adamson_raw'])
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
    adata = adata[:, ~adata.var_names.duplicated()]
    adata.obs = adata.obs[['perturbation', 'is_control', 'perturbation_type']]

    # preprocess  
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=10)

    # - split to inference and evaluation datasets
    ctr_pertb = adata[adata.obs['is_control']].obs['perturbation'].unique()
    non_ctr_pertubs =adata[~adata.obs['is_control']].obs['perturbation'].unique()
    train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.8, random_state=32)
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
    adata_test_bulk = sum_by(adata_test_sc, unique_mapping=True, col='perturbation') 


    # - normalize adata_bulk
    adata_bulk.layers['X_norm'] = adata_bulk.X.copy()

    # - normalize adata_test_sc
    adata_test_sc.layers['X_norm'] = adata_test_sc.X.copy()

    # - normalize evaluation data
    adata_test_bulk.layers['X_norm'] = adata_test_bulk.X.copy()

    # - normalize adata_train_sc
    adata_train_sc.layers['X_norm'] = adata_train_sc.X.copy()

    # - add metadata
    adata_test_sc = add_metadata(adata_test_sc)
    adata_test_bulk = add_metadata(adata_test_bulk)
    adata_train_sc = add_metadata(adata_train_sc)
    adata_bulk = add_metadata(adata_bulk)


    # - save 
    adata_bulk.write(par['adamson_bulk'])
    adata_test_sc.write(par['adamson_test_sc'])
    adata_test_bulk.write(par['adamson_test_bulk'])
    adata_train_sc.write(par['adamson_train_sc'])

def create_test_data(par):
    adata_train_sc = ad.read_h5ad(par['adamson_train_sc'])
    adata_test_bulk = ad.read_h5ad(par['adamson_test_bulk'])
    adata_test_sc = ad.read_h5ad(par['adamson_test_sc'])

    # - subset to 2000 genes
    adata_train_sc_test = adata_train_sc[:, :1000]
    adata_test_bulk_test = adata_test_bulk[:, :1000]
    adata_test_sc_test = adata_test_sc[:, :2000]
    
    # - adata_train_sc: subset to 2 perturbation
    perturbations = adata_train_sc.obs['perturbation'].unique()
    adata_train_sc_test = adata_train_sc[adata_train_sc.obs['perturbation'].isin(perturbations[:2])]
    def sample_per_group(df, n=1000):
        return df.sample(n=min(n, df.shape[0]), replace=False, random_state=42)  # Ensure we don't sample more than available
    sampled_indices = (
        adata_train_sc_test.obs.groupby('perturbation', group_keys=False)
        .apply(sample_per_group)
        .index
    )
    adata_train_sc_test = adata_train_sc_test[sampled_indices]


    # - adata_test_bulk: subset to 2 perturbation
    perturbations = adata_test_bulk.obs['perturbation'].unique()
    adata_test_bulk_test = adata_test_bulk[adata_test_bulk.obs['perturbation'].isin(perturbations[:2])]

    # - adata_test_sc: subset to 2 perturbation
    ctr = adata_test_sc[adata_test_sc.obs['is_control']].obs['perturbation'].unique().tolist()
    perturbations = adata_test_sc.obs['perturbation'].unique().tolist()
    perturbations = list(set(ctr + perturbations[0:1]))  # Ensure unique values
    adata_test_sc_test = adata_test_sc[adata_test_sc.obs['perturbation'].isin(perturbations)]
    def sample_per_group(df, n=1000):
        return df.sample(n=min(n, df.shape[0]), replace=False, random_state=42)  # Ensure we don't sample more than available
    sampled_indices = (
        adata_test_sc_test.obs.groupby('perturbation', group_keys=False)
        .apply(sample_per_group)
        .index
    )
    adata_test_sc_test = adata_test_sc_test[sampled_indices]


    # - add metadata
    adata_train_sc_test = add_metadata(adata_train_sc_test)
    adata_test_bulk_test = add_metadata(adata_test_bulk_test)
    adata_test_sc_test = add_metadata(adata_test_sc_test)


    print(f"adata_train_sc_test: {adata_train_sc_test}")
    print(f"adata_test_bulk_test: {adata_test_bulk_test}")
    print(f"adata_test_sc_test: {adata_test_sc_test}")


    adata_train_sc_test.write(par['adamson_train_sc_test'])
    adata_test_bulk_test.write(par['adamson_test_bulk_test'])
    adata_test_sc_test.write(par['adamson_test_sc_test'])



if __name__ == '__main__':
    main(par)
    create_test_data(par)