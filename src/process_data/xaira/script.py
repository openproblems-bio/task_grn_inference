import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
from sklearn.model_selection import train_test_split 

test_run = True

par = {
    'tf_all': f'resources/grn_benchmark/prior/tf_all.csv',
    'train_share': 0.5,
    
}


raw_dir = '/vol/projects/CIIM/PerturbationDataset/Perturb_seq_dataset_Xiara'

ref_cell_types = ['HEK293T', 'HCT116'] #'HEK293T', 'HCT116' #Human Embryonic Kidney 293T cells, Human Colorectal Carcinoma Cell Line 116


def add_metadata(adata, cell_type):
    adata.uns['dataset_summary'] = f'Whole genome perturbation single cell RNA-seq data on {cell_type}'
    adata.uns['dataset_description'] = 'The rapid expansion of massively parallel sequencing technologies has enabled the development of foundation models to uncover novel biological findings. While these have the potential to significantly accelerate scientific discoveries by creating AI-driven virtual cell models, their progress has been greatly limited by the lack of large-scale high-quality perturbation data, which remains constrained due to scalability bottlenecks and assay variability. Here, we introduce “Fix-Cryopreserve-ScRNAseq” (FiCS) Perturb-seq, an industrialized platform for scalable Perturb-seq data generation. We demonstrate that FiCS Perturb-seq exhibits high sensitivity and low batch effects, effectively capturing perturbation-induced transcriptomic changes and recapitulating known biological pathways and protein complexes. In addition, we release X-Atlas: Orion edition (X-Atlas/Orion), the largest publicly available Perturb-seq atlas. This atlas, generated from two genome-wide FiCS Perturb-seq experiments targeting all human protein-coding genes, comprises eight million cells deeply sequenced to over 16,000 unique molecular identifiers (UMIs) per cell. Furthermore, we show that single guide RNA (sgRNA) abundance can serve as a proxy for gene knockdown (KD) efficacy. Leveraging the deep sequencing and substantial cell numbers per perturbation, we also show that stratification by sgRNA expression can reveal dose-dependent genetic effects. Taken together, we demonstrate that FiCS Perturb-seq is an efficient and scalable platform for high-throughput Perturb-seq screens. Through the release of X-Atlas/Orion, we highlight the potential of FiCS Perturb-seq to address current scalability and variability challenges in data generation, advance foundation model development that incorporates gene-dosage effects, and accelerate biological discoveries.'
    adata.uns['data_reference'] = "@article {Huang2025.06.11.659105, author = {Huang, Ann C and Hsieh, Tsung-Han S and Zhu, Jiang and Michuda, Jackson and Teng, Ashton and Kim, Soohong and Rumsey, Elizabeth M and Lam, Sharon K and Anigbogu, Ikenna and Wright, Philip and Ameen, Mohamed and You, Kwontae and Graves, Christopher J and Kim, Hyunsung John and Litterman, Adam J and Sit, Rene V and Blocker, Alex and Chu, Ci}, title = {X-Atlas/Orion: Genome-wide Perturb-seq Datasets via a Scalable Fix-Cryopreserve Platform for Training Dose-Dependent Biological Foundation Models},elocation-id = {2025.06.11.659105}, year = {2025}, doi = {10.1101/2025.06.11.659105}, publisher = {Cold Spring Harbor Laboratory},URL = {https://www.biorxiv.org/content/early/2025/06/16/2025.06.11.659105}, eprint = {https://www.biorxiv.org/content/early/2025/06/16/2025.06.11.659105.full.pdf}, journal = {bioRxiv}}"
    adata.uns['data_url'] = 'https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1'
    adata.uns['dataset_id'] = f'xaira_{cell_type}'
    adata.uns['dataset_name'] = f'Xaira: {cell_type}'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    return adata
def normalize(adata: ad.AnnData) -> ad.AnnData:
    X_norm = sc.pp.normalize_total(adata, inplace=False)['X']
    X_norm = sc.pp.log1p(X_norm, copy=True)

    adata.layers['X_norm'] = X_norm
    return adata

def split_data(adata: ad.AnnData, train_share):
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    obs = adata.obs
    obs['is_tf'] = obs['perturbation'].isin(tf_all)
    
    unique_perts = obs['perturbation'].unique()
    tf_perturbs = obs[obs['is_tf']]['perturbation'].unique()
    non_tf_perturbs = np.setdiff1d(unique_perts, tf_perturbs)

    # calculate how many TFs and non-TFs to put in train
    n_train_tfs = int(train_share * len(tf_perturbs))
    n_train_non_tfs = int(train_share * len(non_tf_perturbs))

    # sample for train
    np.random.seed(32)
    train_tfs = np.random.choice(tf_perturbs, size=n_train_tfs, replace=False)
    train_non_tfs = np.random.choice(non_tf_perturbs, size=n_train_non_tfs, replace=False)

    # the rest go to test
    test_tfs = np.setdiff1d(tf_perturbs, train_tfs)
    test_non_tfs = np.setdiff1d(non_tf_perturbs, train_non_tfs)

    train_perturbs = np.concatenate([train_tfs, train_non_tfs])
    test_perturbs = np.concatenate([test_tfs, test_non_tfs])

    print(f"Train TFs: {len(train_tfs)}, Train non-TFs: {len(train_non_tfs)}")
    print(f"Test TFs: {len(test_tfs)}, Test non-TFs: {len(test_non_tfs)}")
    print(f"Train total: {len(train_perturbs)}, Test total: {len(test_perturbs)}")

    return train_perturbs, test_perturbs


for ref_cell_type in ref_cell_types:
    print('Reading data for', ref_cell_type, flush=True)
    
    adata = ad.read_h5ad(f'{raw_dir}/{ref_cell_type}_filtered_dual_guide_cells.h5ad', backed='r')
    del adata.obsp 
    del adata.varm
    del adata.uns
    del adata.obsm
    del adata.var 
    del adata.uns

    if test_run: # test
        print('Running in test mode', flush=True)
        test_targets = adata.obs['gene_target'].unique()[:100]
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

        cell_count_t = 20

    # - process single cell data and save it
    adata.obs['is_control'] =  adata.obs['gene_target'] == 'Non-Targeting'
    adata.obs['perturbation_type'] = 'knockdown'
    adata.obs = adata.obs.rename({'gene_target':'perturbation'}, axis=1)[['perturbation', 'is_control', 'perturbation_type', 'sample']]

    # - pseudo bulk 
    if True: 
        print('Creating pseudo bulk data', flush=True)
        from ciim.src.process_dataset.bulkify.helper import bulkify_main
        adata.obs['group'] = np.where(adata.obs['is_control'], adata.obs['sample'], adata.obs['perturbation'])
        adata.obs['group'] = adata.obs['group'].astype('str')
        adata_bulk = bulkify_main(adata, covariates=['group'], cell_count_t=cell_count_t)

    # - train test split
    print('Splitting data...', flush=True)
    train_perturbs, test_perturbs = split_data(adata, par['train_share'])
    
    adata_train_sc = adata[adata.obs['perturbation'].isin(train_perturbs)].copy()
    adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)].copy()
    adata_test_sc = adata[adata.obs['perturbation'].isin(test_perturbs)].copy()


    print('Normalizing...')
    adata_bulk.layers['X_norm'] = adata_bulk.X.copy()
    adata_test_bulk.layers['X_norm'] = adata_test_bulk.X.copy()
   
    adata_train_sc = normalize(adata_train_sc)
    adata_test_sc = normalize(adata_test_sc)
    adata = normalize(adata)
    

    # - add metadata
    adata_bulk = add_metadata(adata_bulk, ref_cell_type)
    adata_test_bulk = add_metadata(adata_test_bulk, ref_cell_type)
    adata = add_metadata(adata, ref_cell_type)
    adata_train_sc = add_metadata(adata_train_sc, ref_cell_type)
    adata_test_sc = add_metadata(adata_test_sc, ref_cell_type)

    print('Saving data...', flush=True)
    adata_train_sc.write_h5ad(f'resources/grn_benchmark/inference_data/xaira_{ref_cell_type}_train_sc.h5ad', compression='gzip')
    adata_test_bulk.write_h5ad(f'resources/grn_benchmark/evaluation_data/xaira_{ref_cell_type}_bulk.h5ad', compression='gzip')
    adata_test_sc.write_h5ad(f'resources/grn_benchmark/evaluation_data/xaira_{ref_cell_type}_sc.h5ad', compression='gzip')

    adata_bulk.write_h5ad(f'resources/extended_data/xaira_{ref_cell_type}_bulk.h5ad', compression='gzip')
    adata.write_h5ad(f'resources/extended_data/xaira_{ref_cell_type}_sc.h5ad', compression='gzip')