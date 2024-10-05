import anndata as ad
import scanpy as sc
## VIASH START
par = {
    # 'multiome_counts': 'resources/datasets_raw/multiome_counts.h5ad',
    # 'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
    # 'multiomics_atac': 'resources/grn-benchmark/multiomics_atac.h5ad'
}
## VIASH END
multiomics_rna = ad.read_h5ad(par['multiomics_rna'])
multiomics_atac = ad.read_h5ad(par['multiomics_atac'])

multiomics_rna_d0 = multiomics_rna[multiomics_rna.obs.donor_id=='donor_0', :]
multiomics_atac_d0 = multiomics_atac[multiomics_atac.obs.donor_id=='donor_0', :]


#  hvgs
var = sc.pp.highly_variable_genes(multiomics_rna_d0, flavor='seurat_v3', n_top_genes=7000, inplace=False)
multiomics_rna_d0.var['highly_variable'] = var.highly_variable

multiomics_rna_d0_hvg = multiomics_rna_d0[:, multiomics_rna_d0.var.highly_variable]

print(multiomics_rna_d0_hvg.shape)
print(multiomics_atac_d0.shape)



multiomics_rna_d0_hvg.write(par['multiomics_rna_d0_hvg'])
multiomics_atac_d0.write(par['multiomics_atac_d0'])