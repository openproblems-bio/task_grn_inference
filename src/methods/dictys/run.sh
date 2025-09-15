
# viash run src/methods/dictys/config.vsh.yaml -- \
# viash run src/methods/dictys/config.vsh.yaml -- ---setup build
# docker run -it --rm -v $(pwd):/workspace -w /workspace ghcr.io/openproblems-bio/task_grn_inference/grn_methods/dictys:dev python src/methods/dictys/script.py \
docker run -it --rm -v $(pwd):/workspace -w /workspace andrewsg/dictys python src/methods/dictys/script.py \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --temp_dir output/temp \
    --prediction output/temp/predictions.h5ad \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv 


