
library(zellkonverter)
options(digits=5, max.print=100)  # Adjust numbers as needed
library(Seurat)
# library(reticulate)


## VIASH START
par <- list(
  perturbation_data_n = "resources_test/grn-benchmark/perturbation_data.h5ad",
  perturbation_data_bc = "resources_test/grn-benchmark/perturbation_data.h5ad"
)
## VIASH END

print(par)
batch_key = 'plate_name'
adata_seurat = readH5AD(par$perturbation_data_n) # raw counts

for (norm_name in c('lognorm', 'pearson')){
  
  
  seurat <- as.Seurat(adata_seurat, counts = "X", data = norm_name)
  batch_list <- SplitObject(seurat, split.by = batch_key)
  anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
  integrated <- IntegrateData(anchors)
  # Extract the integrated expression matrix
  integrated_expr <- GetAssayData(integrated)
  # Make sure the rows and columns are in the same order as the original object
  integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
  # Transpose the matrix to AnnData format
  integrated_expr <- as.matrix(integrated_expr)


  assay(adata_seurat, paste0("seurat_", norm_name)) <- integrated_expr
}

print("Writing adata to h5ad")
writeH5AD(adata_seurat, par$perturbation_data_bc)
