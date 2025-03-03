library(tigress)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
  multiomics_rna = 'resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad',
  tf_all = 'resources/grn_benchmark/prior/tf_all.csv',
  prediction = 'output/tigress_d0_hvg.csv',
  temp_dir = 'output/tigress',
  max_n_links = 50000,
  nsplit = 25
)
## VIASH END

# input expression data
ad <- anndata::read_h5ad(par$multiomics_rna)
X <- as.matrix(ad$X)
gene_names <- colnames(ad)

# Load list of putative TFs
dat <- read.csv(par$tf_all, header = FALSE)
Tf <- intersect(gene_names, dat$V1)

# Run GRN inference method
start.time <- Sys.time()
grn = tigress(
    X, tflist = gene_names, targetlist = gene_names,
    nstepsLARS = 5,
    nsplit = par$nsplit,
    allsteps=FALSE, verb=TRUE, usemulticore=TRUE
)
time.taken <- Sys.time() - start.time

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
