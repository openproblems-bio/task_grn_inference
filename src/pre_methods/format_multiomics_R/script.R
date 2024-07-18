library(anndata)


## VIASH START
par <- list(
  multiomics_rna = "resources_test/grn-benchmark/multiomics_rna.h5ad",
  multiomics_atac = "resources_test/grn-benchmark/multiomics_atac.h5ad",
  temp_dir =  'output/temp_figr/',
  num_workers = 4,
  n_topics = 48,
  prediction= "output/prediction.csv"
)
# meta <- list(
#   functionality_name = "my_method_r"
# )
## VIASH END

cat("Reading input files\n")
input_train <- anndata::read_h5ad(par[["input_train"]])
input_test <- anndata::read_h5ad(par[["input_test"]])

cat("Preprocess data\n")
# ... preprocessing ...

cat("Train model\n")
# ... train model ...

cat("Generate predictions\n")
# ... generate predictions ...

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  obs = list(
    label_pred = obs_label_pred
  ),
  uns = list(
    dataset_id = input_train$uns[["dataset_id"]],
    normalization_id = input_train$uns[["normalization_id"]],
    method_id = meta[["functionality_name"]]
  )
)
output$write_h5ad(par[["output"]], compression = "gzip")