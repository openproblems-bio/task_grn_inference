#!/bin/bash
# viash ns build --parallel
RUN_ID="robust_analy_causal_1" 
resources_dir="resources"
# resources_dir="s3://openproblems-data/resources/grn"
publish_dir="${resources_dir}/results/${RUN_ID}"

reg_type=ridge
subsample=-2
num_workers=10
layer=(scgen_pearson)
metric_ids="[regression_1]"

param_file="./params/${RUN_ID}.yaml"
cat >> $param_file << HERE
param_list:
HERE

# add causal corr
cat >> $param_file << HERE
  - id: corr-causal
    metric_ids: ${metric_ids}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    reg_type: $reg_type
    method_id: baseline_corr_causal
    layer: ${layer}
    subsample: $subsample
    num_workers: $num_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    tf_all: ${resources_dir}/prior/tf_all.csv
    causal: True
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: corr-${1}
    metric_ids: ${metric_ids}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    reg_type: $reg_type
    method_id: baseline_corr-${1}
    layer: ${layer}
    subsample: $subsample
    num_workers: $num_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    tf_all: ${resources_dir}/prior/tf_all.csv
    causal: False
HERE
}


for iter in {1..2}; do  # Loop from 1 to 100 iterations
    append_entry "$iter"   # Execute the append_entry command
done

cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# params_list_file="params/list_${RUN_ID}.yaml"

# param_file="./params/${RUN_ID}.yaml"


# # Loop through grn_names and layers
# layers=("pearson")  # Array containing the layer(s)



# aws s3 sync params/ s3://openproblems-data/resources/grn/params
# # Append the remaining output_state and publish_dir to the YAML file
# cat >> $param_file << HERE
# param_list: "${resources_dir}/${params_list_file}"
# output_state: "state.yaml"
# publish_dir: "$publish_dir"
# HERE

nextflow run . \
  -main-script  target/nextflow/workflows/run_robustness_analysis_causal/main.nf \
  -profile docker \
  -with-trace \
  -c src/common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}
