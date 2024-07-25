import anndata as ad
#!aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/2023-09-14_kaggle_upload/2023-08-31_sc_multiome_expression_atac.h5ad ./resources/raw-data/ --no-sign-request
# mv resources/raw-data/2023-08-31_sc_multiome_expression_atac.h5ad resources/raw-data/multiome.h5ad
par = {
    'multiome_counts': 'resources/raw-data/multiome.h5ad',
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
    'multiomics_atac': 'resources/grn-benchmark/multiomics_atac.h5ad'
}
# Load 
multiomics = ad.read_h5ad(par['multiome_counts'])
multiomics.X = multiomics.layers['counts']
del multiomics.layers
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
multiomics_rna.write(par['multiomics_rna'])
# ATAC
multiomics_atac = multiomics[:,multiomics.var.feature_types=='Peaks']
multiomics_atac.var = multiomics_atac.var[[]]
multiomics_atac.write(par['multiomics_atac'])