library(ennet)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
    "multiomics_rna" = 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
    "tf_all" = 'resources/prior/tf_all.csv',
    "prediction" = 'output/ennet/prediction.csv',
    "temp_dir": 'output/ennet',
    "max_n_links": 50000,
    "M": 100
)
## VIASH END

# input expression data
ad <- anndata::read_h5ad(par$multiomics_rna)
X <- as.matrix(ad$X)
gene_names <- colnames(ad)

# Remove genes with > 90% of zeros
zero_proportion <- colMeans(X == 0)
mask <- (zero_proportion <= 0.9)
X <- X[, mask]
gene_names <- gene_names[mask]

# Remove samples with > 90% of zeros
zero_proportion <- rowMeans(X == 0)
mask <- (zero_proportion <= 0.9)
X <- X[mask,]

# Load list of putative TFs
dat <- read.csv(par$tf_all, header = FALSE)
Tf <- which(gene_names %in% dat$V1)

# Run GRN inference method
K <- matrix(0,nrow(X),ncol(X))
grn <- ennet(E = X, K = K, Tf = Tf, M = par$M)

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
