
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
    'replogle_raw': f'resources/datasets_raw/replogle.h5ad',
    'tf_all': f'resources/grn_benchmark/prior/tf_all.csv',
    'replogle_bulk': f'resources/extended_data/replogle_bulk.h5ad',
    'replogle_test_bulk': f'resources/grn_benchmark/evaluation_data/replogle_bulk.h5ad',
    'replogle_train_bulk': f'resources/grn_benchmark/inference_data/replogle_rna.h5ad',
    'replogle_test_perturbs': f'resources/grn_benchmark/prior/replogle_test_perturbs.csv'
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
    adata.uns['dataset_description'] = 'A central goal of genetics is to define the relationships between genotypes and phenotypes. High-content phenotypic screens such as Perturb-seq (CRISPR-based screens with single-cell RNA-sequencing readouts) enable massively parallel functional genomic mapping but, to date, have been used at limited scales. Here, we perform genome-scale Perturb-seq targeting all expressed genes with CRISPR interference (CRISPRi) across >2.5 million human cells. We use transcriptional phenotypes to predict the function of poorly characterized genes, uncovering new regulators of ribosome biogenesis (including CCDC86, ZNF236, and SPATA5L1), transcription (C7orf26), and mitochondrial respiration (TMEM242). In addition to assigning gene function, single-cell transcriptional phenotypes allow for in-depth dissection of complex cellular phenomena-from RNA processing to differentiation. We leverage this ability to systematically identify genetic drivers and consequences of aneuploidy and to discover an unanticipated layer of stress-specific regulation of the mitochondrial genome. Our information-rich genotype-phenotype map reveals a multidimensional portrait of gene and cellular function.'
    adata.uns['data_reference'] = "@article{replogle2022mapping,\n\ttitle={Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq}, \n\tauthor={Replogle, Joseph M and Saunders, Reuben A and Pogson, Angela N and Hussmann, Jeffrey A and Lenail, Alexander and Guna, Alina and Mascibroda, Lauren and Wagner, Eric J and Adelman, Karen and Lithwick-Yanai, Gila and others},\n\tjournal={Cell},  \n\tvolume={185}, \n\tnumber={14},\n\tpages={2559--2575},\n\tyear={2022},\n\tpublisher={Elsevier}\n\t}"
    adata.uns['data_url'] = 'https://pubmed.ncbi.nlm.nih.gov/35688146/'
    adata.uns['dataset_id'] = 'replogle'
    adata.uns['dataset_name'] = 'Replogle'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'original'
    return adata

if __name__ == '__main__':
    # - get the data
    print('loading...')
    adata = ad.read_h5ad(par['replogle_raw'])
    adata = adata[:, ~adata.var_names.duplicated()]
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
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
    adata_bulk = adata.copy()
    del adata 
    adata_bulk.var.index = adata_bulk.var.index.astype(str)
    adata_bulk.obs = adata_bulk.obs[['perturbation', 'is_control', 'perturbation_type']]

    # preprocess  
    print('preprocessing...')
    sc.pp.filter_cells(adata_bulk, min_genes=100)
    sc.pp.filter_genes(adata_bulk, min_cells=10)


    # - split to inference and evaluation datasets
    if False:
        ctr_pertb = adata_bulk[adata_bulk.obs['is_control']].obs['perturbation'].unique()
        non_ctr_pertubs =adata_bulk[~adata_bulk.obs['is_control']].obs['perturbation'].unique()
        train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.2, random_state=32)
        train_perturbs = np.concatenate([train_perturbs, ctr_pertb]) # add control perturbations to test set for ws_distance
        test_perturbs = np.concatenate([test_perturbs, ctr_pertb]) 
        adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)] 
        adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)] 

    else:
        print('Splitting data...')
        obs = adata_bulk.obs
        obs['is_tf'] = obs['perturbation'].isin(tf_all)
        unique_perts = obs['perturbation'].unique()
        
        tf_all = obs[obs['is_tf']]['perturbation'].unique()
        non_tfs = np.setdiff1d(unique_perts, tf_all)
        train_perturbs, test_perturbs = train_test_split(non_tfs, test_size=.2, random_state=32)

        n_tfs_in_test = int(len(tf_all)/2)

        np.random.seed(32)
        test_tfs = np.random.choice(tf_all, size=n_tfs_in_test, replace=False) 
        train_tfs = np.setdiff1d(tf_all, test_tfs)
        test_perturbs = np.concatenate([test_tfs, test_perturbs])

        train_perturbs = np.concatenate([train_tfs, train_perturbs])

        
        print(f"Test TFs: {n_tfs_in_test}, Train TFs: {len(tf_all) - n_tfs_in_test}")


        adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)]
        adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)]

        print(f"Train : {adata_train_bulk.shape}, {adata_train_bulk.obs['perturbation'].nunique()} \
                Test: {adata_test_bulk.shape}", {adata_test_bulk.obs['perturbation'].nunique()})



    print('filtering...')
    # - filter genes and cells 
    sc.pp.filter_cells(adata_train_bulk, min_genes=100)
    sc.pp.filter_genes(adata_train_bulk, min_cells=10)

    sc.pp.filter_cells(adata_test_bulk, min_genes=100)
    sc.pp.filter_genes(adata_test_bulk, min_cells=10)


    # - normalize 
    print('normalizing...')
    adata_bulk.layers['X_norm'] = adata_bulk.X.copy()
    adata_test_bulk.layers['X_norm'] = adata_test_bulk.X.copy()
    adata_train_bulk.layers['X_norm'] = adata_train_bulk.X.copy()

    # - add metadata
    adata_bulk = add_metadata(adata_bulk)
    adata_test_bulk = add_metadata(adata_test_bulk)
    adata_train_bulk = add_metadata(adata_train_bulk)

    # - save 
    print('saving...')
    adata_bulk.write(par['replogle_bulk'])
    adata_test_bulk.write(par['replogle_test_bulk'])
    adata_train_bulk.write(par['replogle_train_bulk'])

    np.savetxt(par['replogle_test_perturbs'], test_perturbs, fmt="%s",  delimiter=",")
