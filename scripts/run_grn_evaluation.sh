# to run a single GRN evaluation, use the following command:
# bash scripts/run_grn_evaluation.sh --prediction=<prediction file (e.g. prediction.h5ad)> --save_dir=<save dir> --dataset=<dataset (replogle)> --build_images=<true/false (building docker images-only needed one time)> --test_run=<true/false (to use test data)> --no_aws=true


set -e

NO_AWS="true"
RUN_TEST=false
PREDICTION="none"
SAVE_DIR="none"
BUILD_IMAGES=true
reg_type="ridge"
num_workers=4


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
        --no_aws=*)
            NO_AWS="${arg#*=}"
            shift
            ;;
        --reg_type=*)
            reg_type="${arg#*=}"
            shift
            ;;
        --num_workers=*)
            num_workers="${arg#*=}"
            shift
            ;;
        *)
            echo "Unknown argument: $arg"
            exit 1
            ;;
    esac
done

echo "$@"
echo "DATASET: $DATASET"
echo "PREDICTION: $PREDICTION"
echo "RUN_TEST: $RUN_TEST"
echo "BUILD_IMAGES: $BUILD_IMAGES"
echo "NO_AWS: $NO_AWS"
echo "reg_type: $reg_type"

if [ -z "${DATASET:-}" ]; then
    python src/utils/config.py
    source src/utils/config.env
    ALL_DATASETS=(${DATASETS//,/ })
    RUN_ID="all_evaluation"
else
    ALL_DATASETS=("$DATASET")
    RUN_ID="${DATASET}_evaluation"
fi
echo "Datasets to run: ${ALL_DATASETS[*]}"
apply_tf=true
layer='lognorm'
if [ "$RUN_TEST" = "false" ]; then
    resource_folder="resources/"
else
    resource_folder="resources_test/"
fi

if [ "$NO_AWS" = true ]; then
  resources_dir="./${resource_folder}"
else
  resources_dir="s3://openproblems-data/${resource_folder}/grn"
fi

if [ "$SAVE_DIR" != "none" ]; then
  publish_dir="${SAVE_DIR}"
elif [ ${#ALL_DATASETS[@]} -eq 1 ]; then
  publish_dir="${resources_dir}/benchmark/results/${ALL_DATASETS[0]}/"
else
  publish_dir="${resources_dir}/benchmark/results/"
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

if [ "$NO_AWS" = true ]; then
  cat >> "$param_local" << HERE
param_list:
HERE
fi

append_entry() {
  local grn_name="$1"
  local prediction="$2"
  local dataset="$3"

  # Get cell type and metrics from sourced env variables
  cell_type_var="CELLTYPE_${dataset}"
  metrics_var="METRICS_${dataset}"

  cell_type="${!cell_type_var}"
  metric_ids="[${!metrics_var}]"

  echo ${dataset}  ${cell_type} ${metric_ids}

  if [[ "$dataset" =~ ^(norman|nakatake|adamson)$ ]]; then
    layer_='X_norm'
  else
      layer_=$layer
  fi

  cat >> "$param_local" << HERE
  - id: ${grn_name}_${dataset}
    metric_ids: ${metric_ids}
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    prediction: ${prediction}
    num_workers: ${num_workers}
    apply_tf: ${apply_tf}
    reg_type: ${reg_type}
    layer: $layer_
    geneset_hallmark_2020: ${resources_dir}/grn_benchmark/prior/pathways/hallmark_2020.gmt
    geneset_kegg_2021: ${resources_dir}/grn_benchmark/prior/pathways/kegg_2021.gmt
    geneset_reactome_2022: ${resources_dir}/grn_benchmark/prior/pathways/reactome_2022.gmt
    geneset_go_bp_2023: ${resources_dir}/grn_benchmark/prior/pathways/go_bp_2023.gmt
    geneset_bioplanet_2019: ${resources_dir}/grn_benchmark/prior/pathways/bioplanet_2019.gmt
    geneset_wikipathways_2019: ${resources_dir}/grn_benchmark/prior/pathways/wikipathways_2019.gmt
    
HERE
  # Additional fields for specific datasets
  if [[ "$dataset" =~ ^(replogle|xaira_HCT116|xaira_HEK293T|norman)$ ]]; then
    cat >> "$param_local" << HERE
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
HERE
    if [[ "$dataset" =~ ^(replogle|xaira_HCT116|xaira_HEK293T)$ ]]; then
      cat >> "$param_local" << HERE
    evaluation_data_de: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_de.h5ad
HERE
    fi
  fi

  # Add ground truth files for all datasets except nakatake
  if [[ "$dataset" != "nakatake" ]]; then
    cat >> "$param_local" << HERE
    ground_truth_unibind: ${resources_dir}/grn_benchmark/ground_truth/${cell_type}_unibind.csv
    ground_truth_chipatlas: ${resources_dir}/grn_benchmark/ground_truth/${cell_type}_chipatlas.csv
    ground_truth_remap: ${resources_dir}/grn_benchmark/ground_truth/${cell_type}_remap.csv
HERE
  fi
}

if [ "$PREDICTION" != "none" ]; then
  append_entry "single_run" $PREDICTION "${ALL_DATASETS[0]}"
else
  METHODS=(${METHODS//,/ })
  available_methods=()
  for ds in "${ALL_DATASETS[@]}"; do
    grn_models_folder="${resources_dir}/benchmark/results/${ds}/"
    grn_models_folder_local="./resources/benchmark/results/${ds}/"
    for grn_name in "${METHODS[@]}"; do
      prediction_file="${grn_models_folder_local}/${ds}.${grn_name}.${grn_name}.prediction.h5ad"
      if [[ -f "${prediction_file}" ]]; then
        prediction_file=${grn_models_folder}/${ds}.${grn_name}.${grn_name}.prediction.h5ad
        append_entry "$grn_name" $prediction_file "$ds"
        available_methods+=("$grn_name")
      fi
    done
  done
  echo "Available methods:"
  printf '%s\n' "${available_methods[@]}" | sort -u
fi

# Append final fields
if [ "$NO_AWS" = true ]; then
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
    -c common/nextflow_helpers/labels_ci.config
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
    --config scripts/configs/labels_tw.config
fi 