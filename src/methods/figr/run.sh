viash run src/methods/figr/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/ibd_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/ibd_atac.h5ad \
    --prediction resources/results/ibd/ibd.figr.figr.prediction.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --temp_dir output/figr \
    --num_workers 10


