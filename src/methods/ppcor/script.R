library(ppcor)
library(anndata)
library(dplyr)

## VIASH START
par <- list(
  rna = 'resources/grn_benchmark/inference_data/rna.h5ad',
  prediction = 'resources/grn_models/ppcor.h5ad',
  tf_all = 'resources/grn_benchmark/prior/tf_all.csv',
  max_n_links = 50000
)
## VIASH END
args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  if (args[i] == "--rna" && (i+1) <= length(args)) {
    par$rna <- args[i+1]
  } else if (args[i] == "--prediction" && (i+1) <= length(args)) {
    par$prediction <- args[i+1]
  } else if (args[i] == "--layer" && (i+1) <= length(args)) {
    par$layer <- args[i+1]
  }
}



# input expression data
tf_names <- scan(par$tf_all, what = "", sep = "\n")

ad <- anndata::read_h5ad(par$rna)
dataset_id = ad$uns$dataset_id

# Determine which layer to use
if (dataset_id %in% c("nakatake", "norman", "adamson")) {
  layer <- "X_norm"
} else if (!is.null(par$layer)) {
  layer <- par$layer
} else {
  # Default to using X if no layer specified
  layer <- NULL
}

# Get expression data
if (is.null(layer)) {
  inputExpr <- ad$X
} else {
  inputExpr <- ad$layers[[layer]]
  # If the specified layer doesn't exist, throw error
  if (is.null(inputExpr)) {
    stop(paste0("Error in dataset '", dataset_id, "': Layer '", layer, "' not found in AnnData object."))
  }
}

# Additional check to ensure we have data
if (is.null(inputExpr)) {
  stop(paste0("Error in dataset '", dataset_id, "': Could not extract expression data from AnnData object. Both X and specified layer are NULL."))
}

geneNames <- colnames(inputExpr)
colnames(inputExpr) <- c(geneNames)
X <- as.matrix(inputExpr)

# Run GRN inference method
pcorResults = pcor(x = X, method="pearson")

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
print(head(net))
net$weight <- as.character(net$weight)
if (!is.data.frame(net)) {
    stop("Error: 'net' is not a dataframe")
}


output <- AnnData(
  X = matrix(nrow = 0, ncol = 0),
  uns = list(
    method_id = "ppcor",
    dataset_id = dataset_id,
    prediction = net[, c("source", "target", "weight")]
  )
)

print(output)
output$write_h5ad(par$prediction, compression = "gzip")
print("Finished.")
