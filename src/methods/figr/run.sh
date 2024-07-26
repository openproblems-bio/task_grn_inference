viash run src/methods/figr/config.vsh.yaml -- --multiomics_rna resources_test/grn-benchmark/multiomics_r/rna.rds \
 --multiomics_atac resources_test/grn-benchmark/multiomics_r/atac.rds \
 --prediction output/prediction.csv \
 --cell_topic resources_test/grn-benchmark/supp/cell_topic.csv \
 --temp_dir output/figr \
 --num_workers 4




