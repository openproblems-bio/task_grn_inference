viash run src/process_data/nakatake/config.vsh.yaml -- \
    --nakatake_raw resources/datasets_raw/nakatake.h5ad \
    --nakatake_bulk resources/extended_data/nakatake_bulk.h5ad \
    --nakatake_test_bulk resources/grn_benchmark/evaluation_data/nakatake_bulk.h5ad \
    --nakatake_train_bulk resources/grn_benchmark/inference_data/nakatake_rna.h5ad