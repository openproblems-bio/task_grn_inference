

nextflow run target/nextflow/control_methods/pearson_corr/main.nf \
  --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
  --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
  --prediction output/grnboost2.h5ad \
  -profile singularity \
  --publish_dir out