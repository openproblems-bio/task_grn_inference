viash run src/process_data/norman/config.vsh.yaml -- \
    --norman_raw resources/datasets_raw/norman.h5ad \
    --norman_bulk resources/extended_data/norman_bulk.h5ad \
    --norman_test_sc resources/grn_benchmark/evaluation_data/norman_sc.h5ad \
    --norman_test_bulk resources/grn_benchmark/evaluation_data/norman_bulk.h5ad \
    --norman_train_sc resources/grn_benchmark/inference_data/norman_rna.h5ad
