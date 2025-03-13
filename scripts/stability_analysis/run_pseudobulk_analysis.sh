#!/bin/bash
# datasets="norman replogle op nakatake adamson"
datasets=" op "

run_local=false
num_workers=20
metric_ids="[regression_1, regression_2, ws_distance]" #regression_1, regression_2, ws_distance
RUN_ID="pseudobulk_analysis"
reg_type="ridge"
label="pseudobulk_analysis"


grn_names=(
    "pearson_corr"
)

if [ "$run_local" = true ]; then
  resources_dir="./resources/"
else
  resources_dir="s3://openproblems-data/resources/grn"
fi


files_dir="${resources_dir}/grn_benchmark"
publish_dir="${resources_dir}/results/stability_analysis/${RUN_ID}"
echo "Publish dir: $publish_dir"
grn_models_folder="${resources_dir}/grn_models/"
grn_models_folder_local="./resources/grn_models/"


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

granularities=(-1.0 1.0 3.0 5.0 7.0 9.0 11.0 13.0 15.0 17.0 19.0)

append_entry() {
  local grn_name="$1"
  local dataset="$2"
  local granularity="$3"
  
  cat >> "$param_local" << HERE
  - id: ${grn_name}_${dataset}_${granularity}
    metric_ids: ${metric_ids}
    rna: ${files_dir}/inference_data/${dataset}_rna.h5ad
    evaluation_data: ${files_dir}/evaluation_data/${dataset}_bulk.h5ad
    tf_all: ${files_dir}/prior/tf_all.csv
    regulators_consensus: ${files_dir}/prior/regulators_consensus_${dataset}.json
    granularity: ${granularity}
    num_workers: $num_workers

HERE

  # Additional fields for specific datasets
  if [[ "$dataset" == "norman" || "$dataset" == "adamson" || "$dataset" == "replogle" ]]; then
    cat >> "$param_local" << HERE
    evaluation_data_sc: ${files_dir}/evaluation_data/${dataset}_sc.h5ad
    ws_consensus: ${files_dir}/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${files_dir}/prior/ws_distance_background_${dataset}.csv
HERE
  fi
}

# Iterate over datasets and GRN models
for dataset in $datasets; do
  for grn_name in "${grn_names[@]}"; do
    for granularity in "${granularities[@]}"; do
      if [[ -f "${grn_models_folder_local}/${dataset}/${grn_name}.h5ad" ]]; then
        append_entry "$grn_name" "$dataset" "$granularity"
      fi
    done
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
    -main-script  target/nextflow/workflows/run_pseudobulk_analysis/main.nf \
    -profile docker \
    -with-trace \
    -c common/nextflow_helpers/labels_ci.config \
    -params-file ${param_local}
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
    --main-script target/nextflow/workflows/run_pseudobulk_analysis/main.nf \
    --workspace 53907369739130 \
    --compute-env 7gRyww9YNGb0c6BUBtLhDP \
    --params-file ${param_file} \
    --config common/nextflow_helpers/labels_tw.config \
    --labels ${label}
fi 