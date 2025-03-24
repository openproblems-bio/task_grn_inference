viash run src/methods/multi_omics/scglue/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --prediction output/figr_new/prediction.h5ad \
    --temp_dir output/figr_new/ \
    --num_workers 20
