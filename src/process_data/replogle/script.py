import sys
import anndata as ad

import argparse
import numpy as np
from sklearn.model_selection import train_test_split 
import scanpy as sc
import gc
import argparse
import pandas as pd
import os

argument_parser = argparse.ArgumentParser(description='Process Replogle data for GRN inference.')
argument_parser.add_argument('--run_test', action='store_true', help='Run in test mode with a subset of data.')
args = argument_parser.parse_args() 


## VIASH START
par =  {
    'run_test': args.run_test,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv'
}
## VIASH END

meta = {
    'resources_dir': 'src/process_data/'
}
sys.path.append(meta["resources_dir"])

from helper_data import wrapper_large_perturbation_data, split_data_gene_perturbation

def add_metadata(adata):
    adata.uns['dataset_summary'] = 'Single cell RNA-seq data with 231 perturbations (activation) on K562 cells.'
    adata.uns['dataset_description'] = 'A central goal of genetics is to define the relationships between genotypes and phenotypes. High-content phenotypic screens such as Perturb-seq (CRISPR-based screens with single-cell RNA-sequencing readouts) enable massively parallel functional genomic mapping but, to date, have been used at limited scales. Here, we perform genome-scale Perturb-seq targeting all expressed genes with CRISPR interference (CRISPRi) across >2.5 million human cells. We use transcriptional phenotypes to predict the function of poorly characterized genes, uncovering new regulators of ribosome biogenesis (including CCDC86, ZNF236, and SPATA5L1), transcription (C7orf26), and mitochondrial respiration (TMEM242). In addition to assigning gene function, single-cell transcriptional phenotypes allow for in-depth dissection of complex cellular phenomena-from RNA processing to differentiation. We leverage this ability to systematically identify genetic drivers and consequences of aneuploidy and to discover an unanticipated layer of stress-specific regulation of the mitochondrial genome. Our information-rich genotype-phenotype map reveals a multidimensional portrait of gene and cellular function.'
    adata.uns['data_reference'] = "@article{replogle2022mapping,\n\ttitle={Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq}, \n\tauthor={Replogle, Joseph M and Saunders, Reuben A and Pogson, Angela N and Hussmann, Jeffrey A and Lenail, Alexander and Guna, Alina and Mascibroda, Lauren and Wagner, Eric J and Adelman, Karen and Lithwick-Yanai, Gila and others},\n\tjournal={Cell},  \n\tvolume={185}, \n\tnumber={14},\n\tpages={2559--2575},\n\tyear={2022},\n\tpublisher={Elsevier}\n\t}"
    adata.uns['data_url'] = 'https://pubmed.ncbi.nlm.nih.gov/35688146/'
    adata.uns['dataset_id'] = 'replogle'
    adata.uns['dataset_name'] = 'Replogle'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
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

 

def main(par):
    print('Loading single cell data...')
    adata = ad.read_h5ad('resources/datasets_raw/replogle_K562_gwps_raw_singlecell.h5ad', backed='r') 
    if par['run_test']: # to test 
        print('Running in test mode...')
        tf_all = np.loadtxt(par['tf_all'], dtype=str)
        filter = adata.obs['gene'][adata.obs['gene'].isin(tf_all)].unique()[0:100]
        adata = adata[(adata.obs['gene'].isin(filter)) | (adata.obs['gene']=='non-targeting')].to_memory()  
    else:
        print('Running in full mode...')
        adata = adata.to_memory()
        

    # Format the raw data
    print('Formatting raw data...')
    adata = format_raw_data(adata)

    wrapper_large_perturbation_data(adata, split_func=split_data_gene_perturbation,
        covariates=['perturbation'], 
        add_metadata=add_metadata,
        save_name='replogle')
    
    
    print('Data processing complete!')

if __name__ == '__main__':
    main(par)