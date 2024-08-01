
library(scry)
library(zellkonverter)
library(SingleCellExperiment)
options(digits=5, max.print=100)  # Adjust numbers as needed



## VIASH START
par <- list(
  perturbation_data = "resources/grn-benchmark/perturbation_data.h5ad",
  multiomics_rna = "resources/grn-benchmark/multiomics_rna.h5ad",
  hvgs = "resources/grn-benchmark/supp/hvgs.txt",
  n_hvgs = 3000 
)
## VIASH END

print(par)


adata = readH5AD(par$perturbation_data) # raw counts
multiomics_rna <- readH5AD(par$multiomics_rna)

# Extract the gene names from multiomics_rna
multiomics_genes <- rownames(multiomics_rna)

# Subset adata to keep only the genes present in multiomics_rna
adata <- adata[rownames(adata) %in% multiomics_genes, ]

adata_sce = devianceFeatureSelection(adata, assay="X", batch=colData(adata)$plate_name)

binomial_deviance <- rowData(adata_sce)$binomial_deviance

# Sort the indices of binomial deviance in decreasing order and select the top `n_hvgs`
indices <- order(binomial_deviance, decreasing = TRUE)[1:par$n_hvgs]

# Create a mask
mask <- rep(FALSE, length(binomial_deviance))
mask[indices] <- TRUE

# Select the highly variable genes
hvgs_sce <- rownames(adata_sce)[mask]

# Save the highly variable genes to a text file
print(dim(hvgs_sce))

write(hvgs_sce, file = par$hvgs)