
library(zellkonverter)
options(digits=5, max.print=100)  # Adjust numbers as needed
library(Seurat)


## VIASH START
par <- list(
  perturbation_data = "resources_test/grn-benchmark/perturbation_data.h5ad"
)
print(par)
batch_key = 'plate_name'
norm_name = 'lognorm'

adata_seurat = readH5AD(par$perturbation_data) # raw counts



seurat <- as.Seurat(adata_seurat, counts = "X", data = norm_name)
batch_list <- SplitObject(seurat, split.by = batch_key)
anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
integrated <- IntegrateData(anchors)
# Extract the integrated expression matrix
integrated_expr <- GetAssayData(integrated)
# Make sure the rows and columns are in the same order as the original object
integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
# Transpose the matrix to AnnData format
integrated_expr <- t(integrated_expr)