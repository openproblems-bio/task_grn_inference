viash build src/methods/grnboost2/config.vsh.yaml -p docker -o bin/grnboost2 && bin/grnboost2/grnboost2 --multiomics_rna resources/resources_test/grn-benchmark/multiomics_rna.h5ad --prediction output/grnboost2/prediction.csv