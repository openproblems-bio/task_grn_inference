library(ppcor)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
  rna = 'resources/grn-benchmark/rna.h5ad',
  prediction = 'resources/grn_models/ppcor.csv',
  tf_all = 'resources/prior/tf_all.csv',
  max_n_links = 50000
)
## VIASH END
print(par)
# input expression data
tf_names <- scan(par$tf_all, what = "", sep = "\n")

ad <- anndata::read_h5ad(par$rna)

inputExpr <- ad$X
geneNames <- colnames(inputExpr)
colnames(inputExpr) <- c(geneNames)
X <- as.matrix(inputExpr)

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

# Filter to keep only connections where the source is a TF
df_filtered <- df_filtered %>% filter(source %in% tf_names)

# Reset index
df_filtered$index = 1:nrow(df_filtered)

# Keep top links
df_final <- head(df_filtered, par$max_n_links)

# Save results
write.table(df_final, par$prediction, sep = ",", quote = FALSE, row.names = FALSE)

print("Finished.")
