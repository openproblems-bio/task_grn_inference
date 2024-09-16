save_dir="resources/grn_models/donor_0_default"
cell_type_specific="false"

# echo  "negative control"
# viash run src/control_methods/negative_control/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna_0.h5ad \
#   --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
#   --tf_all resources/prior/tf_all.csv \
#   --prediction ${save_dir}/negative_control.csv 


# echo  "baseline pearson"
# viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
#   --tf_all resources/prior/tf_all.csv \
#   --causal false \
#   --corr_method pearson \
#   --cell_type_specific  false \
#   --metacell  false \
#   --impute false \
#   --prediction resources/grn_models/baselines/baseline_pearson.csv 

echo  "baseline pearson causal"
viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna_0.h5ad \
  --tf_all resources/prior/tf_all.csv \
  --causal true \
  --corr_method pearson \
  --cell_type_specific  $cell_type_specific \
  --metacell  false \
  --impute false \
  --prediction ${save_dir}/baseline_pearson_causal.csv 


# echo  "baseline pearson causal metacell"
# viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
#   --tf_all resources/prior/tf_all.csv \
#   --causal true \
#   --corr_method pearson \
#   --cell_type_specific  false \
#   --metacell  true \
#   --impute false \
#   --prediction resources/grn_models/baselines/baseline_pearson_causal_metacell.csv 

# echo  "baseline pearson causal imputation"
# viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
#   --tf_all resources/prior/tf_all.csv \
#   --causal true \
#   --corr_method pearson \
#   --cell_type_specific  false \
#   --metacell  false \
#   --impute true \
#   --prediction resources/grn_models/baselines/baseline_pearson_causal_impute.csv 

# echo  "positive control"
# viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/perturbation_data.h5ad \
#   --tf_all resources/prior/tf_all.csv \
#   --causal true \
#   --corr_method pearson \
#   --cell_type_specific  $cell_type_specific \
#   --metacell  false \
#   --impute false \
#   --prediction ${save_dir}/positive_control.csv 