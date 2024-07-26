viash build src/methods/figr/config.vsh.yaml -p docker -o bin_figr
bin_figr/figr 
 --multiomics_rna resources_test/grn-benchmark/multiomics_r/rna.rds 
 --multiomics_atac resources_test/grn-benchmark/multiomics_r/atac.rds 
 --prediction bin_figr/prediction.csv 
 --cell_topic resources_test/grn-benchmark/supp/cell_topic.csv



