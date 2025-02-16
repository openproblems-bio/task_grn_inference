#!/bin/bash
run_local=false
num_workers=10
metric_ids="[regression_1, regression_2, ws_distance]" #regression_1, regression_2, ws_distance
RUN_ID="scglue_test"
reg_type="ridge"
label=${RUN_ID}

dataset_ids=" op  "
# method_ids="[negative_control, positive_control, pearson_corr, portia, ppcor, scenic, scenicplus , grnboost2]"

method_ids="[scglue]"

echo "Run ID: $RUN_ID"

if [ "$run_local" = true ]; then
  resources_dir="./resources_test/"
else
  resources_dir="s3://openproblems-data/resources_test/grn"
fi

files_dir="${resources_dir}/grn_benchmark"
publish_dir="${resources_dir}/results/${RUN_ID}"
echo "Publish dir: $publish_dir"
grn_models_folder="${resources_dir}/grn_models/"
grn_models_folder_local="./resources/grn_models/"


params_dir="./params"
param_file="${params_dir}/${RUN_ID}.yaml"
param_local="${params_dir}/${RUN_ID}_param_local.yaml"
param_aws="s3://openproblems-data/resources/grn/results/params/${RUN_ID}_param_local.yaml"

# Print GRN names correctly
echo $param_local
echo "GRN models: ${grn_names[@]}"

# Ensure param_file is clean
> "$param_local"
> "$param_file"

if [ "$run_local" = true ]; then
  cat >> "$param_local" << HERE
param_list:
HERE
fi

append_entry() {
  local dataset="$1"
  echo "Dataset: $dataset"
  cat >> "$param_local" << HERE
  - id: ${reg_type}_${grn_name}_${dataset}
    metric_ids: $metric_ids
    method_ids: $method_ids
    rna: ${files_dir}/inference_data/${dataset}_rna.h5ad
    evaluation_data: ${files_dir}/evaluation_data/${dataset}_bulk.h5ad
    tf_all: ${files_dir}/prior/tf_all.csv
    regulators_consensus: ${files_dir}/prior/regulators_consensus_${dataset}.json
    layer: 'X_norm'
    num_workers: $num_workers
    chromsizes: ${files_dir}/prior/chromsizes.csv

HERE
  if [[ "$dataset" == "norman" || "$dataset" == "adamson" || "$dataset" == "replogle" ]]; then
    cat >> "$param_local" << HERE
    evaluation_data_sc: ${files_dir}/evaluation_data/${dataset}_sc.h5ad
    ws_consensus: ${files_dir}/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${files_dir}/prior/ws_distance_background_${dataset}.csv
HERE
  fi
  if [[ "$dataset" == "op" ]]; then
    cat >> "$param_local" << HERE
    atac: ${files_dir}/inference_data/${dataset}_atac.h5ad
HERE
  fi
}

# Iterate over datasets and GRN models
for dataset in $dataset_ids; do
  append_entry "$dataset"
done

# Append final fields
if [ "$run_local" = true ]; then
  cat >> "$param_local" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE
  viash ns build --parallel 
  nextflow run . \
    -main-script  target/nextflow/workflows/run_benchmark/main.nf \
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
    --main-script target/nextflow/workflows/run_benchmark/main.nf \
    --workspace 53907369739130 \
    --compute-env 7gRyww9YNGb0c6BUBtLhDP \
    --params-file ${param_file} \
    --config common/nextflow_helpers/labels_tw.config \
    --labels ${label}
fi

