 
#!/bin/bash

RUN_ID="scores"

echo "Run ID: $RUN_ID"

# resources_dir="./resources/"
resources_dir="s3://openproblems-data/resources/grn"

publish_dir="${resources_dir}/results/${RUN_ID}"
grn_models_folder="${resources_dir}/grn_models/"
grn_models_folder_local="./resources/grn_models/"

reg_type="ridge"
num_workers=10
metric_ids="[regression_1, regression_2, ws_distance]"
datasets="norman adamson"

param_file="./params/${RUN_ID}.yaml"

grn_names=(
    "scglue"
    "scenicplus"
    "celloracle"
    "granie"
    "figr"
    "collectri"
    "grnboost2"
    "ppcor"
    "portia"
    "scenic"
    "positive_control"
    "pearson_corr"
    "negative_control"
    "scprint"
)



# Print GRN names correctly
echo "GRN models: ${grn_names[@]}"

# Ensure param_file is clean
> "$param_file"

# Write YAML header
cat > "$param_file" << HERE
param_list:
HERE

append_entry() {
  local grn_name="$1"
  local dataset="$2"
  
  cat >> "$param_file" << HERE
  - id: ${reg_type}_${grn_name}_${dataset}
    metric_ids: ${metric_ids}
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    reg_type: $reg_type
    method_id: $grn_name
    dataset_id: $dataset
    num_workers: $num_workers
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    prediction: ${grn_models_folder}/${dataset}/${grn_name}.h5ad
    layer: "X_norm"
HERE

  # Additional fields for specific datasets
  if [[ "$dataset" == "norman" || "$dataset" == "adamson" ]]; then
    cat >> "$param_file" << HERE
    evaluation_data_sc: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_sc.h5ad
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
HERE
  fi
}

# Iterate over datasets and GRN models
for dataset in $datasets; do
  for grn_name in "${grn_names[@]}"; do
    if [[ -f "${grn_models_folder_local}/${dataset}/${grn_name}.h5ad" ]]; then
      append_entry "$grn_name" "$dataset"
    fi
  done
done

echo "Parameter file created: $param_file"

# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_grn_inference \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
  --workspace 53907369739130 \
  --compute-env 7gRyww9YNGb0c6BUBtLhDP \
  --params-file ${param_file} \
  --config common/nextflow_helpers/labels_tw.config

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
#   -profile docker \
#   -with-trace \
#   -c common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

# nextflow run openproblems-bio/task_grn_inference -r build/main \
#   -latest \
#   -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
#   -profile singularity \
#   -with-trace \
#   -c common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}


# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_grn_evaluation/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/grn_evaluation_so_ridge.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config


# ./tw launch https://github.com/openproblems-bio/task_grn_inference \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
#   --params-file ${param_file} \
#   --config src/common/nextflow_helpers/labels_tw.config
