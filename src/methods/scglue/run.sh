viash run src/methods/scglue/config.vsh.yaml -- --multiomics_rna resources_test/grn-benchmark/multiomics_rna.h5ad \
    --multiomics_atac resources_test/grn-benchmark/multiomics_atac.h5ad \
    --prediction output/prediction.csv \
    --annotation_file resources_test/grn-benchmark/supp/gencode.v45.annotation.gtf.gz \
    --motif_file resources_test/grn-benchmark/supp/JASPAR2022-hg38.bed.gz \
    --temp_dir output/figr/ \
    --num_workers 2
