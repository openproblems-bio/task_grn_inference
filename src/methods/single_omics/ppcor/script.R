library(ppcor)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
  rna = 'resources/grn-benchmark/rna.h5ad',
  prediction = 'resources/grn_models/ppcor.csv',
  tf_all = 'resources/grn_benchmark/prior/tf_all.csv',
  max_n_links = 50000
)
## VIASH END
args <- commandArgs(trailingOnly = TRUE)

# print(length(args))
# aa
i <- 1
while (i <= length(args)) {
  if (args[i] == "--rna") {
    par$rna <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--prediction") {
    par$prediction <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--tf_all" ) {
    par$tf_all <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--max_n_links") {
    par$max_n_links <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--num_workers") {
    par$num_workers <- as.numeric(args[i + 1])
    i <- i + 2
  } else {
    message("Unknown argument or missing value: ", args[i])
    
  }
}

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
net <- head(df_filtered, par$max_n_links)

# Save results

cat("Output GRN\n")
net$weight <- as.character(net$weight)
output <- AnnData(
  X = matrix(nrow = 0, ncol = 0),
  uns = list(
    method_id = par$method_id,
    dataset_id = par$dataset_id,
    prediction = net[, c("source", "target", "weight")]
  )
)
output$write(par$prediction)

# write.table(df_final, par$prediction, sep = ",", quote = FALSE, row.names = FALSE)

print("Finished.")
