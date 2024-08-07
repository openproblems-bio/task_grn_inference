
library(zellkonverter)
options(digits=5, max.print=100)  # Adjust numbers as needed
library(Seurat)
library(SingleCellExperiment)

# library(reticulate)


## VIASH START
par <- list(
  perturbation_data_n = "resources/grn-benchmark/perturbation_data.h5ad",
  perturbation_data_bc = "resources/grn-benchmark/perturbation_data.h5ad"
)
## VIASH END

print(par)
batch_key = 'plate_name'
adata_seurat = readH5AD(par$perturbation_data_n) # raw counts
print(adata_seurat)
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
  assay_name <- paste0("seurat_", norm_name)
  assays(adata_seurat, withDimnames=FALSE)[[assay_name]] <- integrated_expr
  print(paste0(norm_name, " batch corrected"))
}

print("Writing adata to h5ad")
writeH5AD(adata_seurat, par$perturbation_data_bc)
