library(tigress)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
    "multiomics_rna" = 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
    "prediction" = 'output/tigress/prediction.csv',
    "temp_dir": 'output/tigress',
    "max_n_links": 50000
)
## VIASH END

# input expression data
ad <- anndata::read_h5ad(par$multiomics_rna)
inputExpr <- ad$X

# Run GRN inference method
grn = tigress(inputExpr, allsteps=FALSE, verb=FALSE, usemulticore=TRUE)

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
