


viash run src/methods/granie/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/ibd_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/ibd_atac.h5ad \
    --prediction resources/results/ibd/ibd.granie.granie.prediction.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv
