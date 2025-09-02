viash run src/methods/multi_omics/celloracle/config.vsh.yaml -- --multiomics_rna resources_test/grn-benchmark/multiomics_rna.h5ad \
    --multiomics_atac resources_test/grn-benchmark/multiomics_atac.h5ad \
    --prediction output/prediction.csv \
    --links output/celloracle/links.celloracle.links \
    --base_grn output/celloracle/base_grn.csv \
    --temp_dir output/celloracle/ \
    --num_workers 2
