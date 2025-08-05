#!/bin/bash

# --------------------------
# Dataset-specific availability:
# ws_distance: only for [norman, adamson, replogle]
# scprint: only for [opsca, replogle, norman] (uses different inference data)
# scenicplus, scglue, granie, figr, celloracle: only for [opsca]
# --------------------------

# --- Settings ---
test=false
RUN_ID="nakatake_run"
run_local=false
reg_type="ridge"
num_workers=10
apply_tf_methods=true
apply_skeleton=false

# --- Directories ---
resources_folder=$([ "$test" = true ] && echo "resources_test" || echo "resources")
resources_dir="s3://openproblems-data/${resources_folder}/grn"
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
  local metrics="$2"
  local methods="$3"
  local extra_id="$4"

  cat >> "$param_local" << HERE
  - id="${dataset}${extra_id:+_$extra_id}"
    metric_ids: $metrics
    method_ids: $methods
    rna: ${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad
    rna_all: ${resources_dir}/extended_data/${dataset}_bulk.h5ad
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    layer: 'X_norm'
    num_workers: $num_workers
    apply_tf_methods: $apply_tf_methods
    apply_skeleton: $apply_skeleton
    skeleton: ${resources_dir}/grn_benchmark/prior/skeleton.csv
    reg_type: $reg_type
HERE
  if [ "$extra_id" = "special_case" ]; then # for scprint 
    cat >> "$param_local" << HERE
    rna: ${resources_dir}/grn_benchmark/inference_data/${dataset}_rna_sc_subset.h5ad
HERE
  else  # for rest 
    cat >> "$param_local" << HERE
    rna: ${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad
HERE
  fi

  # Extra fields for certain datasets
  if [[ "$dataset" =~ ^(norman|adamson|replogle)$ ]]; then
    cat >> "$param_local" << HERE
    evaluation_data_sc: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_sc.h5ad
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
HERE
  fi

  if [[ "$dataset" =~ ^(op|opsca)$ ]]; then
    cat >> "$param_local" << HERE
    atac: ${resources_dir}/grn_benchmark/inference_data/${dataset}_atac.h5ad
HERE
  fi
}

# --------- COMBINATIONS TO ADD ----------

# append_entry "op" "[regression_1,regression_2, ws_distance]" "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, ppcor, scenic, scprint, grnboost,
#                                                                         scenicplus, scglue, granie, figr, celloracle]" 
# append_entry "norman"  "[regression_1,regression_2, ws_distance]" "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, ppcor, scenic, scprint, grnboost]"
# append_entry "adamson"  "[regression_1,regression_2, ws_distance]" "[pearson_corr, negative_control, positive_control, 
#                                                                         portia, ppcor, scenic, grnboost]"
append_entry "nakatake"  "[regression_1,regression_2]" "[pearson_corr, negative_control, positive_control, 
                                                                        portia, scenic, grnboost]"
# append_entry "replogle" "[regression_1, regression_2, ws_distance]" "[pearson_corr, negative_control, positive_control, portia, ppcor, scenic, grnboost]"
# append_entry "replogle" "[regression_1, regression_2, ws_distance]" "[scprint]" "special_case"                                                
# --- Final configuration ---
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

  aws s3 cp $param_local $param_aws

  # tw launch https://github.com/openproblems-bio/task_grn_inference \
  #   --revision build/main \
  #   --pull-latest \
  #   --main-script target/nextflow/workflows/run_benchmark/main.nf \
  #   --workspace 53907369739130 \
  #   --compute-env 7gRyww9YNGb0c6BUBtLhDP \
  #   --params-file ${param_file} \
  #   --config common/nextflow_helpers/labels_tw.config \
  #   --labels ${RUN_ID}
  tw launch https://github.com/openproblems-bio/task_grn_inference \
    --revision build/main \
    --pull-latest \
    --main-script target/nextflow/workflows/run_benchmark/main.nf \
    --workspace 53907369739130 \
    --compute-env 6TJs9kM1T7ot4DbUY2huLF \
    --params-file ${param_file} \
    --labels ${RUN_ID} \
    --config scripts/labels_tw.config
fi

#on demand 6TJs9kM1T7ot4DbUY2huLF   
#7gRyww9YNGb0c6BUBtLhDP