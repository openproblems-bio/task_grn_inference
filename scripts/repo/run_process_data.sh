
RUN_ID="process_data"
echo $RUN_ID
raw_resources_dir="./resources/datasets_raw/"
publish_dir="./output/process_data"

param_file="./params/${RUN_ID}.yaml"

echo $param_file
# Start writing to the YAML file
cat > $param_file << HERE
param_list:
  - id: all_datasets
    op_multiome: ${raw_resources_dir}/op_multiome_sc_counts.h5ad
    op_perturbation_raw: ${raw_resources_dir}/op_perturbation_sc_counts.h5ad
publish_dir: "$publish_dir"
output_state: "state.yaml"
HERE

viash ns build --parallel 

nextflow run . \
  -main-script  target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  -with-trace \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}