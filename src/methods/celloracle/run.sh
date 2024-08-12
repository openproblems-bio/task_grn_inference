viash run src/methods/celloracle/config.vsh.yaml -- --multiomics_rna resources_test/grn-benchmark/multiomics_rna.h5ad \
    --multiomics_atac resources_test/grn-benchmark/multiomics_atac.h5ad \
    --prediction output/prediction.csv \
    --temp_dir output/celloracle/ \
    --num_workers 2
