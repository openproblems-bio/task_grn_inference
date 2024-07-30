import anndata as ad

par = {
    # 'multiome_counts': 'resources/datasets_raw/multiome_counts.h5ad',
    # 'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
    # 'multiomics_atac': 'resources/grn-benchmark/multiomics_atac.h5ad'
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

# ATAC
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

multiomics_rna.write(par['multiomics_rna'])
multiomics_atac.write(par['multiomics_atac'])