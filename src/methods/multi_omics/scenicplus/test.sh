viash run src/methods/multi_omics/scenicplus/config.vsh.yaml -- --multiomics_atac resources_test/grn-benchmark/multiomics_atac.h5ad \
    --multiomics_rna resources_test/grn-benchmark/multiomics_rna.h5ad \
    --temp_dir output/scenicplus \
    --prediction output/scenicplus/prediction.csv \
    --cell_topic output/scenicplus/cell_topic.csv \
    --scplus_mdata output/scenicplus/scplus_mdata.h5mu
