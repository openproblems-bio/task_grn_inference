import anndata as ad
import scanpy as sc
import numpy as np
## VIASH START
par = {
    'multiome_counts': 'resources/datasets_raw/op_multiome_sc_counts.h5ad',
    'multiomics_rna': 'resources/grn_benchmark/inference_data/op_rna.h5ad',
    'multiomics_atac': 'resources/grn_benchmark/inference_data/op_atac.h5ad'
}
## VIASH END

if __name__ == '__main__':

    # Load 
    multiomics = ad.read_h5ad(par['multiome_counts'])
    multiomics.X = multiomics.layers['counts']
    del multiomics.layers
    multiomics.layers['counts'] = multiomics.X.copy()

    multiomics.var.index.name='location'
    multiomics.obs.index.name='obs_id'

    # map the cell types
    cell_types_o = multiomics.obs.cell_type.unique()
    T_cell_types = ['T regulatory cells', 'T cells CD8+', 'T cells CD4+']
    cell_type_map = {cell_type: 'T cells' if cell_type in T_cell_types else cell_type for cell_type in cell_types_o}
    multiomics.obs['cell_type'] = multiomics.obs['cell_type'].map(cell_type_map)
    # RNA
    multiomics_rna = multiomics[:,multiomics.var.feature_types=='Gene Expression']
    multiomics_rna.var = multiomics_rna.var[['gene_ids', 'interval']]

    # def high_coverage(adata):
    #     threshold = 0.1
    #     mask = adata.X!=0
    #     mask_obs = (np.sum(mask, axis=1).A.flatten()/mask.shape[1])>threshold
    #     mask_var = (np.sum(mask, axis=0).A.flatten()/mask.shape[0])>threshold
    #     adata.obs['high_coverage'] = mask_obs
    #     adata.var['high_coverage'] = mask_var
    # high_coverage(multiomics_rna)
    #------ ATAC
    multiomics_atac = multiomics[:,multiomics.var.feature_types=='Peaks']
    multiomics_atac.var = multiomics_atac.var[[]]

    # Find common cells (observations) in both RNA and ATAC datasets
    common_obs = multiomics_rna.obs_names.intersection(multiomics_atac.obs_names)
    # Subset the RNA and ATAC data to keep only the common cells
    multiomics_rna = multiomics_rna[common_obs, :]
    multiomics_atac = multiomics_atac[common_obs, :]

    print(multiomics_rna)
    print(multiomics_atac)

    # change donor names
    unique_donors = multiomics_rna.obs.donor_id.unique()
    donor_map = {donor_id: f'donor_{i}' for i, donor_id in enumerate(unique_donors)}
    multiomics_rna.obs['donor_id'] = multiomics_rna.obs['donor_id'].map(donor_map)
    multiomics_atac.obs['donor_id'] = multiomics_atac.obs['donor_id'].map(donor_map)

    # normalize rna 
    X_norm = sc.pp.normalize_total(multiomics_rna, inplace=False)['X']
    multiomics_rna.layers['X_norm'] = sc.pp.log1p(X_norm, copy=True)

    multiomics_rna.uns['dataset_id'] = 'op'
    multiomics_atac.uns['dataset_id'] = 'op'

    multiomics_rna.uns['dataset_summary'] = 'RNA-seq data from the OPSCA dataset'
    multiomics_atac.uns['dataset_summary'] = 'ATAC-seq data from the OPSCA dataset'

    multiomics_rna.uns['dataset_organism'] = 'human'
    multiomics_atac.uns['dataset_organism'] = 'human'

    multiomics_rna.uns['normalization_id'] = 'sla'
    multiomics_atac.uns['normalization_id'] = 'sla'


    multiomics_rna.write(par['multiomics_rna'])
    multiomics_atac.write(par['multiomics_atac'])