viash build src/methods/single_omics/genie3/config.novsh.yaml -p docker -o bin/genie3 && bin/genie3/genie3 --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad --tfs resources/prior/tf_all.csv --prediction output/genie3/prediction.csv