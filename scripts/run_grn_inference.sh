#!/bin/bash

set -e 

# --- Settings ---
RUN_TEST=false
num_workers=20
apply_tf_methods=true
layer='lognorm'
NO_AWS=false
# Parse arguments
for arg in "$@"; do
    case $arg in
        --dataset=*)
            DATASET="${arg#*=}"
            shift
            ;;
        --test_run=*)
            RUN_TEST="${arg#*=}"
            shift
            ;;
        --no_aws=*)
            NO_AWS="${arg#*=}"
            shift
            ;;
        *)
            echo "Unknown argument: $arg"
            exit 1
            ;;
    esac
done
if [ -z "${DATASET:-}" ]; then
    python src/utils/config.py
    source src/utils/config.env
    ALL_DATASETS=(${DATASETS//,/ })
    RUN_ID="all_inference"
else
    ALL_DATASETS=("$DATASET")
    RUN_ID="${DATASET}_inference"
fi
echo "Datasets to run: ${ALL_DATASETS[*]}"

# --- Directories ---
resources_folder=$([ "$RUN_TEST" = true ] && echo "resources_test" || echo "resources")
if [ "$NO_AWS" = true ]; then
  resources_dir="./${resources_folder}/"
else
  resources_dir="s3://openproblems-data/${resources_folder}/grn"
fi

if [ ${#ALL_DATASETS[@]} -eq 1 ]; then
    publish_dir="${resources_dir}/benchmark/results/${ALL_DATASETS[0]}"
else
    publish_dir="${resources_dir}/benchmark/results/"
fi


params_dir="./params"
param_file="${params_dir}/${RUN_ID}.yaml"
param_local="${params_dir}/${RUN_ID}_param_local.yaml"
param_aws="s3://openproblems-data/resources/grn/results/params/${RUN_ID}_param_local.yaml"

echo "Resources folder: $resources_folder"
echo "Run ID: $RUN_ID"
echo "Publish dir: $publish_dir"
echo "Local param file: $param_local"

# --- Prepare param files ---
# Ensure param_file is clean
> "$param_local"
> "$param_file"

if [ "$NO_AWS" = true ]; then
  cat >> "$param_local" << HERE
param_list:
HERE
fi


# --- Function to append dataset block ---
append_entry() {
  local dataset="$1"
  local methods="$2"
  local use_train_sc=false

  # check if third argument is non-empty (or truthy)
  if [ -n "$3" ]; then
      use_train_sc=true
  fi

  if [[ "$dataset" =~ ^(norman|nakatake|adamson)$ ]]; then
    layer_='X_norm'
  else
    layer_="$layer"
  fi
  group_id="${dataset}"
  if [ "$use_train_sc" = true ]; then
    rna_file="${resources_dir}/extended_data/${dataset}_train_sc.h5ad"
  else
    rna_file="${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad"
    
  fi

  cat >> "$param_local" << HERE
  - id: ${group_id}
    method_ids: $methods
    rna: $rna_file
    rna_all: ${resources_dir}/extended_data/${dataset}_rna_all.h5ad
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    layer: $layer_
    num_workers: $num_workers
    apply_tf_methods: $apply_tf_methods
    temp_dir: $dataset
HERE

  if [[ "$dataset" =~ ^(op|opsca|ibd_cd|ibd_uc)$ ]]; then
    cat >> "$param_local" << HERE
    atac: ${resources_dir}/grn_benchmark/inference_data/${dataset}_atac.h5ad
HERE
  fi
}

for DATASET in "${ALL_DATASETS[@]}"; do
  if [[ "$DATASET" =~ ^(replogle|parsebioscience|xaira_HEK293T|xaira_HCT116)$ ]]; then
    methods="[pearson_corr, negative_control, positive_control, grnboost, portia, scenic, geneformer, scgpt, spearman_corr]"
    # methods="[pearson_corr, negative_control, positive_control, portia, geneformer, scgpt, spearman_corr]"
    append_entry "$DATASET" "$methods"
    append_entry "$DATASET" "[scprint]" "true"
    echo $methods
  elif [ "$DATASET" = "op" ] || [ "$DATASET" = "ibd_cd" ] || [ "$DATASET" = "ibd_uc" ]; then
    methods="[pearson_corr, spearman_corr, negative_control, positive_control, grnboost, portia, scenic, scprint, geneformer, scgpt, figr, scenicplus, celloracle, granie, scglue]"
    # methods="[celloracle, scglue]"
    append_entry "$DATASET" "$methods"
    echo $methods
  else
    methods="[pearson_corr, negative_control, positive_control, grnboost, portia, scenic, scprint, geneformer, scgpt, spearman_corr]"
    append_entry "$DATASET" "$methods"
    echo $methods
  fi
done


if [ "$RUN_TEST" = true ]; then
  labels_config="scripts/configs/labels_tw_test.config"
else
  labels_config="scripts/configs/labels_tw.config"
fi
echo "Labels config: $labels_config"
# --- Final configuration ---
if [ "$NO_AWS" = true ]; then
  cat >> "$param_local" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

  viash ns build --parallel 
  nextflow run . \
    -main-script  target/nextflow/workflows/run_grn_inference/main.nf \
    -profile docker \
    -with-trace \
    -c $labels_config \
    -params-file ${param_local}

else
  cat >> "$param_file" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
param_list: "$param_aws"
HERE

  aws s3 cp $param_local $param_aws
  # echo "Launching task_grn_inference on hpc compute..."
  # tw launch https://github.com/openproblems-bio/task_grn_inference \
  #     --revision build/main \
  #     --pull-latest \
  #     --main-script target/nextflow/workflows/run_grn_inference/main.nf \
  #     --workspace 209741690280743 \
  #     --params-file ${param_file} \
  #     --labels ${RUN_ID} \
  #     --config $labels_config
  echo "Launching task_grn_inference on aws compute..."
  tw launch https://github.com/openproblems-bio/task_grn_inference \
    --revision build/main \
    --pull-latest \
    --main-script target/nextflow/workflows/run_grn_inference/main.nf \
    --workspace 53907369739130 \
    --compute-env 6TJs9kM1T7ot4DbUY2huLF \
    --params-file ${param_file} \
    --labels ${RUN_ID} \
    --config $labels_config
fi

#on demand 6TJs9kM1T7ot4DbUY2huLF   
#7gRyww9YNGb0c6BUBtLhDP