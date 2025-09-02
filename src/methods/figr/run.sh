viash run src/methods/multi_omics/figr/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --prediction output/prediction.h5ad \
    --cell_topic resources_test/grn_benchmark/prior/cell_topic.csv \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
    --temp_dir output/figr \
    --num_workers 4

