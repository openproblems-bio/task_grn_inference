import sys
import anndata as ad

import argparse
import numpy as np
from sklearn.model_selection import train_test_split 
import scanpy as sc
import gc

## VIASH START
par =  {
    'replogle_gwps': 'resources/datasets_raw/replogle_K562_gwps_raw_singlecell.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'replogle_gwps_test_sc': 'resources/grn_benchmark/evaluation_data/replogle_sc.h5ad',
    'replogle_gwps_train_sc': 'resources/extended_data/replogle_train_sc.h5ad',
    'replogle_gwps_train_sc_subset': 'resources/grn_benchmark/inference_data/replogle_rna_sc_subset.h5ad',
    'replogle_test_perturbs': f'resources/grn_benchmark/prior/replogle_test_perturbs.csv',
}
## VIASH END
if False:
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--input', type=str, required=True)
    argparser.add_argument('--tf_all', type=str, required=True)
    argparser.add_argument('--adata_test_sc', type=str, required=True)
    argparser.add_argument('--adata_train_sc', type=str, required=True)
    argparser.add_argument('--adata_train_sc_subset', type=str, required=True)

    args = argparser.parse_args()
    par = vars(args)

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        'resources_dir': 'src/utils/'
    }
    sys.path.append(meta["resources_dir"])

from util import sum_by, fetch_gene_info

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
    test_perturbs = np.loadtxt(par['replogle_test_perturbs'], dtype=str)
    test_tfs = np.intersect1d(tf_all, test_perturbs)
    train_pertubs = np.setdiff1d(unique_perts, test_perturbs)

    ctrs =  ['non-targeting']
    test_tfs = np.concatenate([test_tfs, ctrs])
    
    print(f"Test TFs: {len(test_tfs)}")

    adata_test = adata[adata.obs['perturbation'].isin(test_tfs)]
    adata_train = adata[adata.obs['perturbation'].isin(train_pertubs)] # all perturbs

    adata_train_sc_subset = adata[adata.obs['perturbation'].isin(train_pertubs[:500])]

    print(f"Train : {adata_train.shape}, Test: {adata_test.shape}, Train_short: {adata_train_sc_subset.shape}")

    return adata_train, adata_test, adata_train_sc_subset



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
    adata = ad.read_h5ad(par['replogle_gwps'], backed='r') 
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
        adata_bulk.write(par['replogle_gwps_bulk'])
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
        adata_train_bulk.write(par['replogle_gwps_train_bulk'])
        
        
        del adata_train_bulk
        gc.collect()

    adata_train_sc = add_metadata(adata_train_sc)
    adata_train_sc.write(par['replogle_gwps_train_sc'])

    adata_train_sc_subset = adata_train_sc_subset.to_memory()
    adata_train_sc_subset = add_metadata(adata_train_sc_subset)
    adata_train_sc_subset.write(par['replogle_gwps_train_sc_subset'])
    del adata_train_sc, adata_train_sc_subset
    gc.collect()

    # - process test data
    print('Process test data...')
    adata_test_sc = adata_test_sc.to_memory()
    adata_test_sc = normalize(adata_test_sc)
    adata_test_sc = add_metadata(adata_test_sc)
    adata_test_sc.uns['normalization_id'] = 'sla'

    adata_test_sc.write(par['replogle_gwps_test_sc'])

    if False:
        adata_test_bulk = sum_by(adata_test_sc, col='perturbation', unique_mapping=True)
        adata_test_bulk.uns['dataset_id'] = 'replogle'
        adata_test_bulk = normalize_pearson(adata_test_bulk)
        adata_test_bulk.write(par['replogle_gwps_test_bulk'])


if __name__ == '__main__':
    main(par)

