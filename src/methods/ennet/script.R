library(ennet)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
    "multiomics_rna" = 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
    "prediction" = 'output/ennet/prediction.csv',
    "temp_dir": 'output/ennet',
    "max_n_links": 50000
)
## VIASH END

# input expression data
ad <- anndata::read_h5ad(par$multiomics_rna)
X <- ad$X

# Run GRN inference method
K <- matrix(0,nrow(X),ncol(X))
Tf <- 1:ncol(X)
grn = ennet(E = X, K = K, Tf = Tf)

# Re-format output
df <- as.data.frame(as.table(grn))
colnames(df) <- c("source", "target", "weight")
df <- df[df$weight != 0,]
df <- df[order(-df$weight),]

# Add index as extra column
df <- cbind(index = 1:nrow(df), df)

# Keep top links
df <- head(df, par$max_n_links)

# Save results
write.table(df, par$prediction, sep = ",", quote = FALSE, row.names = FALSE)

print("Finished.")
