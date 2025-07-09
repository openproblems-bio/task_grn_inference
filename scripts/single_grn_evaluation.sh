# bash src/metrics/all_metrics/run.sh resources/grn_models/norman/grnboost2.csv norman

prediction=${1}
dataset=${2}
test_run=False
# Parse --test_run flag
for arg in "$@"; do
  if [ "$arg" == "--test_run" ]; then
    test_run=true
  fi
done


run_on_seqera=false
# Parse --run_on_seqera flag
for arg in "$@"; do
  if [ "$arg" == "--run_on_seqera" ]; then
    run_on_seqera=true
  fi
done

if [ -z "$prediction" ]; then
  echo "Error: prediction argument is missing (pass it as the first argument)"
  exit 1
fi

if [ -z "$dataset" ]; then
  echo "Error: dataset name is missing (pass it as the second argument)"
  exit 1
fi

mkdir -p output

RUN_ID="single_evaluation"

publish_dir="./output/${RUN_ID}"

num_workers=1

param_file="./output/${RUN_ID}.yaml"

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

if [ "$test_run" = true ]; then
  echo "Running on test data..."
  resources_dir="./resources_test/" # change this to resources/ for the actual scores
  num_workers=1
else
  echo "Running on actual data..."
  resources_dir="./resources/" 
fi

if [ "$run_on_seqera" = true ]; then
  prediction_aws="s3://openproblems-data/resources/grn/output/prediction_test.h5ad"
  aws s3 cp $prediction $prediction_aws
  prediction=$prediction_aws  

  resources_dir="s3://openproblems-data/resources/grn"
  publish_dir="s3://openproblems-data/resources/grn/output/${RUN_ID}"
fi

echo "Publish dir: $publish_dir"

append_entry() {
  cat >> $param_file << HERE
  - id: ${reg_type}
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    reg_type: "ridge"
    num_workers: $num_workers
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    prediction: $prediction
HERE
  # Additional fields for specific datasets
  if [[ "$dataset" == "norman" || "$dataset" == "adamson" || "$dataset" == "replogle" ]]; then
    cat >> "$param_file" << HERE
    evaluation_data_sc: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_sc.h5ad
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
HERE
  fi
}

append_entry 

# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

echo "Generated parameter file: $param_file"
# build the images: this can be skipped after the first run


if [ "$run_on_seqera" = true ]; then
  echo "Running on Seqera platform..."

  tw launch https://github.com/openproblems-bio/task_grn_inference \
      --revision build/main \
      --pull-latest \
      --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
      --workspace 209741690280743 \
      --params-file ${param_file} \
      --labels ${RUN_ID} \
      --config scripts/hpc_settings.config
  # tw launch https://github.com/openproblems-bio/task_grn_inference \
  #   --revision build/main \
  #   --pull-latest \
  #   --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
  #   --workspace 53907369739130 \
  #   --compute-env 7gRyww9YNGb0c6BUBtLhDP \
  #   --params-file ${param_file} \
  #   --labels ${RUN_ID} \
  #   --config common/nextflow_helpers/labels_tw.config
else
  echo "Running locally..."
  viash ns build --parallel --setup build -s src/metrics/

  viash ns build --parallel 
  nextflow run . \
  -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
  -profile docker \
  -with-trace \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}
fi
