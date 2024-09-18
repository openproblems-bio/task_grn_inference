viash run src/control_methods/baseline_corr/config.vsh.yaml -- \
  --prediction output/baseline_corr.csv \
  --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
  --tf_all resources/prior/tf_all.csv \
  --causal true
