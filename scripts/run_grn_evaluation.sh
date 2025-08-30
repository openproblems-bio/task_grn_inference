# to run a single GRN evaluation, use the following command:
# bash scripts/run_grn_evaluation.sh --prediction=<prediction file (e.g. prediction.h5ad)> --save_dir=<save dir> --dataset=<dataset (replogle)> --build_images=<true/false (building docker images-only needed one time)> --test_run=<true/false (to use test data)> --run_local=true>



#!/bin/bash
set -e

RUN_LOCAL="true"
RUN_TEST=false
PREDICTION="none"
SAVE_DIR="none"
BUILD_IMAGES=true


# Parse arguments
for arg in "$@"; do
    case $arg in
        --dataset=*)
            DATASET="${arg#*=}"
            shift
            ;;
        --prediction=*)
            PREDICTION="${arg#*=}"
            shift
            ;;
        --test_run=*)
            RUN_TEST="${arg#*=}"
            shift
            ;;
        --save_dir=*)
            SAVE_DIR="${arg#*=}"
            shift
            ;;
        --build_images=*)
            BUILD_IMAGES="${arg#*=}"
            shift
            ;;
        --run_local=*)
            RUN_LOCAL="${arg#*=}"
            shift
            ;;
        *)
            echo "Unknown argument: $arg"
            exit 1
            ;;
    esac
done

if [ -z "${DATASET:-}" ]; then
    echo "Error: DATASET must be provided. Use --dataset=<dataset_name>."
    exit 1
fi

num_workers=10
metric_ids="[regression_1, regression_2, ws_distance]" #regression_1, regression_2, ws_distance
RUN_ID="${DATASET}_evaluation"

models_folder="${DATASET}/"
reg_type="ridge"
apply_skeleton=false
apply_tf=true
layer='lognorm'
if [ "$RUN_TEST" = "false" ]; then
    resource_folder="resources/"
else
    resource_folder="resources_test/"
fi

if [ "$RUN_LOCAL" = true ]; then
  resources_dir="./${resource_folder}"
else
  resources_dir="s3://openproblems-data/${resource_folder}/grn"
fi

if [ "$SAVE_DIR" != "none" ]; then
  publish_dir="${SAVE_DIR}"
else
  publish_dir="${resources_dir}/results/${models_folder}"
fi

mkdir -p "$publish_dir"
echo "Publish dir: $publish_dir"


params_dir="./params"
mkdir -p "$params_dir"
param_file="${params_dir}/${RUN_ID}.yaml"
param_local="${params_dir}/${RUN_ID}_param_local.yaml"
param_aws="s3://openproblems-data/resources/grn/results/params/${RUN_ID}_param_local.yaml"

# Ensure param_file is clean
> "$param_local"
> "$param_file"

if [ "$RUN_LOCAL" = true ]; then
  cat >> "$param_local" << HERE
param_list:
HERE
fi

# Write YAML header

append_entry() {
  local grn_name="$1"
  local prediction="$2"
  local dataset="$3"
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
    prediction: ${prediction}
    skeleton: ${resources_dir}/grn_benchmark/prior/skeleton.csv
    apply_skeleton: ${apply_skeleton}
    apply_tf: ${apply_tf}
    layer: $layer_

HERE

  # Additional fields for specific datasets
  if [[ "$dataset" =~ ^(norman|replogle|adamson|xaira_HCT116|xaira_HEK293T)$ ]]; then
    cat >> "$param_local" << HERE
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
HERE
  fi
}

if [ "$PREDICTION" != "none" ]; then
  append_entry "single_run" $PREDICTION "$DATASET"
else
  grn_names=(
      "positive_control"
      "pearson_corr"
      "negative_control"
      "scglue"
      "scenicplus"
      "celloracle"
      "granie"
      "figr"
      "grnboost"
      "ppcor"
      "portia"
      "scenic"
      "scprint"
  )
  grn_models_folder="${resources_dir}/results/${models_folder}/"
  grn_models_folder_local="./resources/results/${models_folder}/" # just to control the hetergenity of the models for different datasets

  # Iterate over GRN models
  available_methods=()
  for grn_name in "${grn_names[@]}"; do
    prediction_file="${grn_models_folder_local}/${DATASET}.${grn_name}.${grn_name}.prediction.h5ad"
    if [[ -f "${prediction_file}" ]]; then
      prediction_file=${grn_models_folder}/${DATASET}.${grn_name}.${grn_name}.prediction.h5ad
      append_entry "$grn_name" $prediction_file "$DATASET"
      available_methods+=("$grn_name")
    fi
  done
  echo "Available methods:"
  printf '%s\n' "${available_methods[@]}" | sort -u

fi


# Append final fields
if [ "$RUN_LOCAL" = true ]; then
  cat >> "$param_local" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE
  if [ "$BUILD_IMAGES" = true ]; then
    echo "Building Docker images..."
    viash ns build --parallel --setup build -s src/metrics/
  else
    viash ns build --parallel 
  fi
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
  # echo "Launching task_grn_inference on aws compute..."
  tw launch https://github.com/openproblems-bio/task_grn_inference \
    --revision build/main \
    --pull-latest \
    --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
    --workspace 53907369739130 \
    --params-file ${param_file} \
    --labels ${RUN_ID} \
    --config common/nextflow_helpers/labels_tw.config
fi 