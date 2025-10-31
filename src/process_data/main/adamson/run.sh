viash run src/process_data/adamson/config.vsh.yaml -- \
    --adamson_raw resources/datasets_raw/adamson.h5ad \
    --adamson_bulk resources/extended_data/adamson_bulk.h5ad \
    --adamson_test_sc resources/grn_benchmark/evaluation_data/adamson_sc.h5ad \
    --adamson_test_bulk resources/grn_benchmark/evaluation_data/adamson_bulk.h5ad \
    --adamson_train_sc resources/grn_benchmark/inference_data/adamson_rna.h5ad \
    --adamson_train_sc_test resources_test/grn_benchmark/inference_data/adamson_rna.h5ad \
    --adamson_test_bulk_test resources_test/grn_benchmark/evaluation_data/adamson_bulk.h5ad \
    --adamson_test_sc_test resources_test/grn_benchmark/evaluation_data/adamson_sc.h5ad