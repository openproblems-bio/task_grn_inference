#!/bin/bash

# --------------------------
# Dataset-specific availability:
# ws_distance: only for [norman, adamson, replogle]
# scprint: only for [opsca, replogle, norman] (uses different inference data)
# scenicplus, scglue, granie, figr, celloracle: only for [opsca]
# --------------------------
set -e 
# --- Settings ---
test=false
RUN_ID="replogle_run"
run_local=false
num_workers=10
apply_tf_methods=true
layer='lognorm'

# --- Directories ---
resources_folder=$([ "$test" = true ] && echo "resources_test" || echo "resources")
if [ "$run_local" = true ]; then
  resources_dir="./${resources_folder}/"
else
  resources_dir="s3://openproblems-data/${resources_folder}/grn"
fi

publish_dir="${resources_dir}/results/${RUN_ID}"
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

if [ "$run_local" = true ]; then
  cat >> "$param_local" << HERE
param_list:
HERE
fi

# --- Function to append dataset block ---
append_entry() {
  local dataset="$1"
  local methods="$2"

  if [[ "$dataset" =~ ^(norman|nakatake|adamson)$ ]]; then
    layer_='X_norm'
  else
      layer_=$layer
  fi
  
  cat >> "$param_local" << HERE
  - id: ${dataset}
    method_ids: $methods
    rna: ${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad
    rna_all: ${resources_dir}/extended_data/${dataset}_bulk.h5ad
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    layer: $layer_
    num_workers: $num_workers
    apply_tf_methods: $apply_tf_methods
HERE

  if [[ "$dataset" =~ ^(op|opsca)$ ]]; then
    cat >> "$param_local" << HERE
    atac: ${resources_dir}/grn_benchmark/inference_data/${dataset}_atac.h5ad
HERE
  fi
}

# --------- COMBINATIONS TO ADD ----------

# append_entry "op" "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, ppcor, scenic, scprint, grnboost,
#                                                                         scenicplus, scglue, granie, figr, celloracle]" 
# append_entry "norman"  "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, ppcor, scenic, scprint, grnboost]"
# append_entry "adamson" "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, ppcor, scenic, grnboost]"
# append_entry "nakatake"  "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, scenic, grnboost]"
# append_entry "replogle" "[pearson_corr, negative_control, positive_control, portia, ppcor, scenic, grnboost, scprint]"

append_entry "replogle" "[pearson_corr, negative_control, positive_control, scprint]" 
# append_entry "xaira_HEK293T" "[pearson_corr, negative_control, positive_control]" 
# append_entry "xaira_HCT116" "[pearson_corr, negative_control, positive_control]"
# append_entry "parsebioscience" "[pearson_corr, negative_control, positive_control]"


# --- Final configuration ---
if [ "$run_local" = true ]; then
  cat >> "$param_local" << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

  viash ns build --parallel 
  nextflow run . \
    -main-script  target/nextflow/workflows/run_grn_inference/main.nf \
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

  aws s3 cp $param_local $param_aws

  
  
  # tw launch https://github.com/openproblems-bio/task_grn_inference \
  #     --revision build/main \
  #     --pull-latest \
  #     --main-script target/nextflow/workflows/run_grn_inference/main.nf \
  #     --workspace 209741690280743 \
  #     --params-file ${param_file} \
  #     --labels ${RUN_ID} \
  #     --config scripts/hpc_settings.config

  tw launch https://github.com/openproblems-bio/task_grn_inference \
    --revision build/main \
    --pull-latest \
    --main-script target/nextflow/workflows/run_grn_inference/main.nf \
    --workspace 53907369739130 \
    --compute-env 6TJs9kM1T7ot4DbUY2huLF \
    --params-file ${param_file} \
    --labels ${RUN_ID} \
    --config scripts/labels_tw.config
fi

#on demand 6TJs9kM1T7ot4DbUY2huLF   
#7gRyww9YNGb0c6BUBtLhDP