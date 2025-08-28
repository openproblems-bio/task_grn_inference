
# viash run src/methods/multi_omics/dictys/config.vsh.yaml -- \
#     --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
#     --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
#     --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
#     --prediction output/d_test/prediction.h5ad \
#     --temp_dir output/d_test/ \
#     --num_workers 10


# python src/methods/multi_omics/dictys/greta/extract_data.py \
#     --pre_path output/dictys/mudata.h5mu \
#     --exp_path output/dictys/exp.csv \
#     --pks_path output/dictys/peak.csv 

bash src/methods/multi_omics/dictys/tfb.sh \
    --input_pre output/dictys/mudata.h5mu \
    --output_d output/dictys \
    --input_frags output/dictys/donor_0.tsv.gz \
    --input_motif output/d_test/data/motif_file.motif \
    --input_genome output/d_test/data/genome \
    --threads 10

