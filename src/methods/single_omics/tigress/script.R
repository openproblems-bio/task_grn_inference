library(tigress)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
    "multiomics_rna" = 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
    "tfs" = 'resources/prior/tf_all.csv',
    "prediction" = 'output/tigress/prediction.csv',
    "temp_dir": 'output/tigress',
    "max_n_links": 50000
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
colnames(X) <- gene_names

# Remove samples with > 90% of zeros
zero_proportion <- rowMeans(X == 0)
mask <- (zero_proportion <= 0.9)
X <- X[mask,]

# Load list of putative TFs
dat <- read.csv(par$tfs, header = FALSE)
Tf <- intersect(gene_names, dat$V1)

# Run GRN inference method
grn = tigress(X, tflist = Tf, targetlist = gene_names, allsteps=FALSE, verb=FALSE, usemulticore=TRUE)

# Re-format output
df <- as.data.frame(as.table(grn))
colnames(df) <- c("source", "target", "weight")
df <- df[df$weight != 0,]
df <- df[order(-df$weight),]

# Add index as extra column
df <- cbind(index = 1:nrow(df), df)

# Keep top links
df <- head(df, par$max_n_links)

print(df)

# Save results
write.table(df, par$prediction, sep = ",", quote = FALSE, row.names = FALSE)

print("Finished.")
