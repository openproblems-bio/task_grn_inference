import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys
import argparse

argument_parser = argparse.ArgumentParser(description='Process Replogle data for GRN inference.')
argument_parser.add_argument('--run_test', action='store_true', help='Run in test mode with a subset of data.')
args = argument_parser.parse_args() 

meta = {
    'resources_dir': 'src/process_data/'
}
sys.path.append(meta["resources_dir"])

from helper_data import wrapper_large_perturbation_data, split_data_gene_perturbation

par = {
    'run_test': args.run_test    
}


raw_dir = '/vol/projects/CIIM/PerturbationDataset/Perturb_seq_dataset_Xiara'

ref_cell_types = ['HEK293T', 'HCT116'] #'HEK293T', 'HCT116' #Human Embryonic Kidney 293T cells, Human Colorectal Carcinoma Cell Line 116

    
def add_metadata(cell_type):     
    def inner(adata):
        adata.uns['dataset_summary'] = f'Whole genome perturbation single cell RNA-seq data on {cell_type}'
        adata.uns['dataset_description'] = 'The rapid expansion of massively parallel sequencing technologies has enabled the development of foundation models to uncover novel biological findings. While these have the potential to significantly accelerate scientific discoveries by creating AI-driven virtual cell models, their progress has been greatly limited by the lack of large-scale high-quality perturbation data, which remains constrained due to scalability bottlenecks and assay variability. Here, we introduce “Fix-Cryopreserve-ScRNAseq” (FiCS) Perturb-seq, an industrialized platform for scalable Perturb-seq data generation. We demonstrate that FiCS Perturb-seq exhibits high sensitivity and low batch effects, effectively capturing perturbation-induced transcriptomic changes and recapitulating known biological pathways and protein complexes. In addition, we release X-Atlas: Orion edition (X-Atlas/Orion), the largest publicly available Perturb-seq atlas. This atlas, generated from two genome-wide FiCS Perturb-seq experiments targeting all human protein-coding genes, comprises eight million cells deeply sequenced to over 16,000 unique molecular identifiers (UMIs) per cell. Furthermore, we show that single guide RNA (sgRNA) abundance can serve as a proxy for gene knockdown (KD) efficacy. Leveraging the deep sequencing and substantial cell numbers per perturbation, we also show that stratification by sgRNA expression can reveal dose-dependent genetic effects. Taken together, we demonstrate that FiCS Perturb-seq is an efficient and scalable platform for high-throughput Perturb-seq screens. Through the release of X-Atlas/Orion, we highlight the potential of FiCS Perturb-seq to address current scalability and variability challenges in data generation, advance foundation model development that incorporates gene-dosage effects, and accelerate biological discoveries.'
        adata.uns['data_reference'] = "@article {Huang2025.06.11.659105, author = {Huang, Ann C and Hsieh, Tsung-Han S and Zhu, Jiang and Michuda, Jackson and Teng, Ashton and Kim, Soohong and Rumsey, Elizabeth M and Lam, Sharon K and Anigbogu, Ikenna and Wright, Philip and Ameen, Mohamed and You, Kwontae and Graves, Christopher J and Kim, Hyunsung John and Litterman, Adam J and Sit, Rene V and Blocker, Alex and Chu, Ci}, title = {X-Atlas/Orion: Genome-wide Perturb-seq Datasets via a Scalable Fix-Cryopreserve Platform for Training Dose-Dependent Biological Foundation Models},elocation-id = {2025.06.11.659105}, year = {2025}, doi = {10.1101/2025.06.11.659105}, publisher = {Cold Spring Harbor Laboratory},URL = {https://www.biorxiv.org/content/early/2025/06/16/2025.06.11.659105}, eprint = {https://www.biorxiv.org/content/early/2025/06/16/2025.06.11.659105.full.pdf}, journal = {bioRxiv}}"
        adata.uns['data_url'] = 'https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1'
        adata.uns['dataset_id'] = f'xaira_{cell_type}'
        adata.uns['dataset_name'] = f'Xaira: {cell_type}'
        adata.uns['dataset_organism'] = 'human'
        adata.uns['normalization_id'] = 'lognorm'
        return adata
    return inner

for ref_cell_type in ref_cell_types:
    print('Reading data for', ref_cell_type, flush=True)
    
    adata = ad.read_h5ad(f'{raw_dir}/{ref_cell_type}_filtered_dual_guide_cells.h5ad', backed='r')
    

    if par['run_test']: # test
        print('Running in test mode', flush=True)
        test_targets = adata.obs['gene_target'].unique()[:20]
        # Initialize a boolean mask of all False
        mask = pd.Series(False, index=adata.obs_names)

        for gene in test_targets:
            gene_cells = adata.obs.index[adata.obs['gene_target'] == gene]
            n_sample = min(10, len(gene_cells))
            sampled_cells = np.random.choice(gene_cells, size=n_sample, replace=False)
            mask.loc[sampled_cells] = True

        adata = adata[mask].to_memory()
        cell_count_t = 1
    else:
        # - QC
        print('Sending to memory', flush=True)
        adata = adata.to_memory()
        print('Running QC', flush=True)
        adata = adata[(adata.obs['n_genes_by_counts']>10) & (adata.obs['n_genes_by_counts']<5000) & (adata.obs['pct_counts_mt']<10)] 
        n_batches = adata.obs['sample'].nunique()
        min_cells = 10*n_batches
        sc.pp.filter_genes(adata, min_cells=min_cells)
        cell_count_t = 10

    # - process single cell data and save it
    adata.obs['is_control'] =  adata.obs['gene_target'] == 'Non-Targeting'
    adata.obs['perturbation_type'] = 'knockdown'
    adata.obs = adata.obs.rename({'gene_target':'perturbation'}, axis=1)[['perturbation', 'is_control', 'perturbation_type', 'sample']]

    # - 
    wrapper_large_perturbation_data(adata, split_func=split_data_gene_perturbation, save_name=f'xaira_{ref_cell_type}', 
                    add_metadata=add_metadata(ref_cell_type), covariates=['perturbation'], cell_count_t=cell_count_t)
    