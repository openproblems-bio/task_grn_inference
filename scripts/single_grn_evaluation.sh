# bash src/metrics/all_metrics/run.sh resources/grn_models/norman/grnboost2.csv norman

prediction=${1}
dataset=${2}

if [ -z "$prediction" ]; then
  echo "Error: prediction argument is missing"
  exit 1
fi

if [ -z "$dataset" ]; then
  echo "Error: dataset argument is missing"
  exit 1
fi

mkdir -p output

RUN_ID="test_run"
echo $RUN_ID
resources_dir="./resources/"
publish_dir="./output/${RUN_ID}"

grn_models_folder="${resources_dir}/grn_models/"

reg_type="ridge"
num_workers=4
metric_ids="[regression_1, regression_2, ws_distance]"

param_file="./output/${RUN_ID}.yaml"

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: ${reg_type}
    metric_ids: ${metric_ids}
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    reg_type: $reg_type
    num_workers: $num_workers
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    prediction: $prediction
    layer: "X_norm"
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


# Loop through grn_names and layers

append_entry 


# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# build the images: this can be skipped after the first run
viash ns build --parallel --setup build -s src/metrics/

viash ns build --parallel 

nextflow run . \
  -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
  -profile docker \
  -with-trace \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}