library(ppcor)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
    "multiomics_rna" = 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
    "prediction" = 'output/ppcor/prediction.csv',
    "temp_dir": 'output/ppcor',
    "max_n_links": 50000
)
## VIASH END

# input expression data
ad <- anndata::read_h5ad(par$multiomics_rna)
inputExpr <- ad$X
geneNames <- colnames(inputExpr)
colnames(inputExpr) <- c(geneNames)
X <- as.matrix(inputExpr)

# Keep genes with less than 10% of zeros
mask <- (apply(X, 2, function(x) mean(x != 0)) >= 0.1)
X <- X[, mask]
geneNames <- geneNames[mask]

# Keep samples with less than 10% of zeros
mask <- (apply(X, 1, function(x) mean(x != 0)) >= 0.1)
X <- X[mask,]

# Run GRN inference method
pcorResults = pcor(x = X, method = "pearson")

# Save results
df = data.frame(
    index = 1:length(pcorResults$estimate),
    source = geneNames[c(row(pcorResults$estimate))],
    target = geneNames[c(col(pcorResults$estimate))],
    weight = c(pcorResults$estimate)
)
df <- df[order(df$weight, decreasing=TRUE),]

# Remove self-interactions
df_filtered <- df[df$source != df$target,]

# Reset index
df_filtered$index = 1:nrow(df_filtered)

# Keep top links
df_final <- head(df_filtered, par$max_n_links)

# Save results
write.table(df_final, par$prediction, sep = ",", quote = FALSE, row.names = FALSE)

print("Finished.")
