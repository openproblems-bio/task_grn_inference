
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
    'nakatake_raw': f'resources/datasets_raw/nakatake.h5ad',
    'nakatake_bulk': f'resources/extended_data/nakatake_bulk.h5ad',
    'nakatake_test_bulk': f'resources/grn_benchmark/evaluation_data/nakatake_bulk.h5ad',
    'nakatake_train_bulk': f'resources/grn_benchmark/inference_data/nakatake_rna.h5ad'
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
    adata.uns['dataset_summary'] = 'RNA-seq data with 463 perturbations (overexpression) on SEES3 cells.'
    adata.uns['dataset_description'] = 'Transcription factors (TFs) play a pivotal role in determining cell states, yet our understanding of the causative relationship between TFs and cell states is limited. Here, we systematically examine the state changes of human pluripotent embryonic stem cells (hESCs) by the large-scale manipulation of single TFs. We establish 2,135 hESC lines, representing three clones each of 714 doxycycline (Dox)-inducible genes including 481 TFs, and obtain 26,998 microscopic cell images and 2,174 transcriptome datasets-RNA sequencing (RNA-seq) or microarrays-48 h after the presence or absence of Dox. Interestingly, the expression of essentially all the genes, including genes located in heterochromatin regions, are perturbed by these TFs. TFs are also characterized by their ability to induce differentiation of hESCs into specific cell lineages. These analyses help to provide a way of classifying TFs and identifying specific sets of TFs for directing hESC differentiation into desired cell types.'
    adata.uns['data_reference'] = "@article{nakatake2020generation, \n\ttitle={Generation and profiling of 2,135 human ESC lines for the systematic analyses of cell states perturbed by inducing single transcription factors}, \n\tauthor={Nakatake, Yuhki and Ko, Shigeru BH and Sharov, Alexei A and Wakabayashi, Shunichi and Murakami, Miyako and Sakota, Miki and Chikazawa, Nana and Ookura, Chiaki and Sato, Saeko and Ito, Noriko and others}, \n\tjournal={Cell Reports},\n\tvolume={31}, \n\tnumber={7}, \n\tyear={2020}, \n\tpublisher={Elsevier}}"
    adata.uns['data_url'] = 'https://pubmed.ncbi.nlm.nih.gov/32433964/'
    adata.uns['dataset_id'] = 'nakatake'
    adata.uns['dataset_name'] = 'Nakatake'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'original'
    return adata

def main(par):
    # - get the data
    adata = ad.read_h5ad(par['nakatake_raw'])
    adata = adata[:, ~adata.var_names.duplicated()]
    adata = adata.copy()
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
    nakatake_bulk = adata.copy()
    del adata 
    nakatake_bulk.var.index = nakatake_bulk.var.index.astype(str)
    nakatake_bulk.obs = nakatake_bulk.obs[['perturbation', 'is_control', 'perturbation_type']]

    # preprocess  
    sc.pp.filter_cells(nakatake_bulk, min_genes=100)
    sc.pp.filter_genes(nakatake_bulk, min_cells=10)


    # - split to inference and evaluation datasets
    ctr_pertb = nakatake_bulk[nakatake_bulk.obs['is_control']].obs['perturbation'].unique()
    non_ctr_pertubs =nakatake_bulk[~nakatake_bulk.obs['is_control']].obs['perturbation'].unique()
    train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.5, random_state=32)
    train_perturbs = np.concatenate([train_perturbs, ctr_pertb]) # add control perturbations to test set for ws_distance
    test_perturbs = np.concatenate([test_perturbs, ctr_pertb]) 

    nakatake_train_bulk = nakatake_bulk[nakatake_bulk.obs['perturbation'].isin(train_perturbs)] 
    nakatake_test_bulk = nakatake_bulk[nakatake_bulk.obs['perturbation'].isin(test_perturbs)] 


    # - filter genes and cells 
    sc.pp.filter_cells(nakatake_train_bulk, min_genes=100)
    sc.pp.filter_genes(nakatake_train_bulk, min_cells=10)

    sc.pp.filter_cells(nakatake_test_bulk, min_genes=100)
    sc.pp.filter_genes(nakatake_test_bulk, min_cells=10)


    # - normalize 
    nakatake_bulk.layers['X_norm'] = nakatake_bulk.X.copy()
    nakatake_test_bulk.layers['X_norm'] = nakatake_test_bulk.X.copy()
    nakatake_train_bulk.layers['X_norm'] = nakatake_train_bulk.X.copy()

    # - add metadata
    nakatake_bulk = add_metadata(nakatake_bulk)
    nakatake_test_bulk = add_metadata(nakatake_test_bulk)
    nakatake_train_bulk = add_metadata(nakatake_train_bulk)
    
    # - save 
    print('Saving...')
    nakatake_bulk.write(par['nakatake_bulk'])
    nakatake_test_bulk.write(par['nakatake_test_bulk'])
    nakatake_train_bulk.write(par['nakatake_train_bulk'])

if __name__ == '__main__':
    main(par)