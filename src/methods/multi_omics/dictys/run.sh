
viash run src/methods/multi_omics/dictys/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
    --prediction output/d_test/prediction.h5ad \
    --temp_dir output/d_test/ \
    --num_workers 10

