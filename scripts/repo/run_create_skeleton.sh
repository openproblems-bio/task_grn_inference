viash run src/stability_analysis/skeleton/config.vsh.yaml -- \
    --atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --annotation_file resources/supp_data/gencode.v45.annotation.gtf.gz \
    --temp_dir output/skeleton \
    --extend_range 150000 \
    --flank_length 1000 \
    --skeleton resources/grn_benchmark/prior/skeleton.csv \
    --motif_dataset_encode resources/supp_data/databases/scglue/ENCODE-TF-ChIP-hg38.bed.gz \
    --motif_dataset_jaspar resources/supp_data/databases/scglue/JASPAR2022-hg38.bed.gz