viash run src/methods/single_omics/grnboost2/config.vsh.yaml -- --multiomics_rna resources_test/grn-benchmark/multiomics_rna.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --cell_type_specific false \
    --prediction output/grnboost2/prediction.csv