viash run src/process_data/replogle/config.vsh.yaml -- \
    --replogle_raw resources/datasets_raw/replogle.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --replogle_bulk resources/extended_data/replogle_bulk.h5ad \
    --replogle_test_bulk resources/grn_benchmark/evaluation_data/replogle_bulk.h5ad \
    --replogle_train_bulk resources/grn_benchmark/inference_data/replogle_rna.h5ad \
    --replogle_test_perturbs resources/grn_benchmark/prior/replogle_test_perturbs.csv