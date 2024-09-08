#!/bin/bash
# viash ns build --parallel
RUN_ID="robust_analy_causal" 
# resources_dir="resources"
resources_dir="s3://openproblems-data/resources/grn"

publish_dir="${resources_dir}/results/${RUN_ID}"

reg_type=ridge
subsample=-2
max_workers=10

param_file="./params/${RUN_ID}.yaml"
# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: corr-${1}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    reg_type: $reg_type
    method_id: corr-${1}
    layer: ${2}
    subsample: $subsample
    max_workers: $max_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    tf_all: ${resources_dir}/prior/tf_all.csv
HERE
}
# Loop through grn_names and layers
layers=("pearson")  # Array containing the layer(s)

for layer in "${layers[@]}"; do  # Iterate over each layer in the array
    for iter in {1..10}; do  # Loop from 1 to 100 iterations
        append_entry "$iter" "$layer"  # Execute the append_entry command
    done
done


# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_robustness_analysis_causal/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

