import anndata as ad
import scanpy as sc
import numpy as np
## VIASH START
par = {
    'op_multiome': 'resources/datasets_raw/op_multiome_sc_counts.h5ad',
    'op_rna': 'resources/grn_benchmark/inference_data/op_rna.h5ad',
    'op_atac': 'resources/grn_benchmark/inference_data/op_atac.h5ad'
}
## VIASH END

if __name__ == '__main__':

    # Load 
    multiomics = ad.read_h5ad(par['op_multiome'])
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
    rna = multiomics[:,multiomics.var.feature_types=='Gene Expression']
    rna.var = rna.var[['gene_ids', 'interval']]

    # def high_coverage(adata):
    #     threshold = 0.1
    #     mask = adata.X!=0
    #     mask_obs = (np.sum(mask, axis=1).A.flatten()/mask.shape[1])>threshold
    #     mask_var = (np.sum(mask, axis=0).A.flatten()/mask.shape[0])>threshold
    #     adata.obs['high_coverage'] = mask_obs
    #     adata.var['high_coverage'] = mask_var
    # high_coverage(rna)
    #------ ATAC
    atac = multiomics[:,multiomics.var.feature_types=='Peaks']
    atac.var = atac.var[[]]

    # Find common cells (observations) in both RNA and ATAC datasets
    common_obs = rna.obs_names.intersection(atac.obs_names)
    # Subset the RNA and ATAC data to keep only the common cells
    rna = rna[common_obs, :]
    atac = atac[common_obs, :]

    print(rna)
    print(atac)

    # change donor names
    unique_donors = rna.obs.donor_id.unique()
    donor_map = {donor_id: f'donor_{i}' for i, donor_id in enumerate(unique_donors)}
    rna.obs['donor_id'] = rna.obs['donor_id'].map(donor_map)
    atac.obs['donor_id'] = atac.obs['donor_id'].map(donor_map)

    # normalize rna 
    X_norm = sc.pp.normalize_total(rna, inplace=False)['X']
    rna.layers['X_norm'] = sc.pp.log1p(X_norm, copy=True)

    rna.uns['dataset_id'] = 'op'
    atac.uns['dataset_id'] = 'op'

    rna.uns['dataset_name'] = 'OPSCA'
    atac.uns['dataset_id'] = 'OPSCA'

    rna.uns['dataset_summary'] = 'RNA-seq data from the OPSCA dataset'
    atac.uns['dataset_summary'] = 'ATAC-seq data from the OPSCA dataset'

    rna.uns['dataset_organism'] = 'human'
    atac.uns['dataset_organism'] = 'human'

    rna.uns['normalization_id'] = 'sla'
    atac.uns['normalization_id'] = 'sla'

    # - needed for some R packages
    annotation_peak = atac.var.reset_index().location.str.split(':', expand=True)
    atac.var['seqname'] = annotation_peak[0].values
    atac.var['ranges'] = annotation_peak[1].values
    atac.var['strand'] = '+'


    rna.write(par['op_rna'])
    atac.write(par['op_atac'])