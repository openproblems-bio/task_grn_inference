#!/bin/bash
DATASET="${1:-replogle}" # Default dataset if not provided
run_local="${2:-false}" 

# datasets="norman replogle op nakatake adamson"
datasets="$DATASET" #xaira_HCT116 xaira_HEK293T parsebioscience replogle

num_workers=10
metric_ids="[regression_1, regression_2, ws_distance]" #regression_1, regression_2, ws_distance
RUN_ID=$DATASET
models_folder="${DATASET}/"
reg_type="ridge"
apply_skeleton=false
apply_tf=true
layer='lognorm'

grn_names=(
    "positive_control"
    "pearson_corr"
    "negative_control"
    "scglue"
    "scenicplus"
    "celloracle"
    "granie"
    "figr"
    "grnboost2"
    "ppcor"
    "portia"
    "scenic"
    "scprint"
)

if [ "$run_local" = true ]; then
  resources_dir="./resources/"
else
  resources_dir="s3://openproblems-data/resources/grn"
fi


publish_dir="${resources_dir}/results/${models_folder}"
echo "Publish dir: $publish_dir"
grn_models_folder="${resources_dir}/results/${models_folder}/"
grn_models_folder_local="./resources/results/${models_folder}/" # just to control the hetergenity of the models for different datasets


params_dir="./params"
param_file="${params_dir}/${RUN_ID}.yaml"
param_local="${params_dir}/${RUN_ID}_param_local.yaml"
param_aws="s3://openproblems-data/resources/grn/results/params/${RUN_ID}_param_local.yaml"

# Print GRN names correctly
echo "GRN models: ${grn_names[@]}"

# Ensure param_file is clean
> "$param_local"
> "$param_file"

if [ "$run_local" = true ]; then
  cat >> "$param_local" << HERE
param_list:
HERE
fi

# Write YAML header

append_entry() {
  local grn_name="$1"
  local dataset="$2"
  if [[ "$dataset" =~ ^(norman|nakatake|adamson)$ ]]; then
    layer_='X_norm'
  else
      layer_=$layer
  fi
  cat >> "$param_local" << HERE
  - id: ${reg_type}_${grn_name}_${dataset}
    metric_ids: ${metric_ids}
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    prediction: ${grn_models_folder}/${dataset}.${grn_name}.${grn_name}.prediction.h5ad
    skeleton: ${resources_dir}/grn_benchmark/prior/skeleton.csv
    apply_skeleton: ${apply_skeleton}
    apply_tf: ${apply_tf}
    layer: $layer_

HERE

  # Additional fields for specific datasets
  if [[ "$dataset" =~ ^(norman|replogle|adamson|xaira_HCT116|xaira_HEK293T)$ ]]; then
    cat >> "$param_local" << HERE
    evaluation_data_sc: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_sc.h5ad
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
HERE
  fi
}

# Iterate over datasets and GRN models
for dataset in $datasets; do
  for grn_name in "${grn_names[@]}"; do
    prediction_file="${grn_models_folder_local}/${dataset}.${grn_name}.${grn_name}.prediction.h5ad"
    if [[ -f "${prediction_file}" ]]; then
      append_entry "$grn_name" "$dataset"
    else
      echo "File not found: ${prediction_file}"
    fi
  done
done


# Append final fields
if [ "$run_local" = true ]; then
  cat >> "$param_local" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE
  viash ns build --parallel 
  echo "Parameter file created: $param_local"
  nextflow run . \
    -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
    -profile docker \
    -with-trace \
    -params-file ${param_local} \
    # -c common/nextflow_helpers/labels_ci.config
else
  cat >> "$param_file" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
param_list: "$param_aws"
HERE
  echo "Parameter file created: $param_file"

  aws s3 cp $param_local $param_aws

  tw launch https://github.com/openproblems-bio/task_grn_inference \
    --revision build/main \
    --pull-latest \
    --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
    --workspace 53907369739130 \
    --params-file ${param_file} \
    --labels ${RUN_ID} \
    --config common/nextflow_helpers/labels_tw.config
    # --config scripts/hpc_settings.config \
    
    # --config common/nextflow_helpers/labels_tw.config
fi 