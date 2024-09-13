viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
  --tf_all resources/prior/tf_all.csv \
  --causal true \
  --prediction output/baseline_causal.csv 

viash run src/metrics/regression_1/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --prediction output/baseline_causal.csv \
    --score output/score_causal.h5ad

viash run src/control_methods/baseline_corr/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
  --tf_all resources/prior/tf_all.csv \
  --causal false \
  --prediction output/baseline_noncausal.csv 

viash run src/metrics/regression_1/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --prediction output/baseline_noncausal.csv \
    --score output/score_noncausal.h5ad