#!/bin/bash

degrees=(0 10 20 50 100)
noise_type="$1" #"net", "weight", "sign"
echo $noise_type

RUN_ID="robust_analy_$1" 
# resources_dir="resources"
resources_dir="s3://openproblems-data/resources/grn"

publish_dir="${resources_dir}/results/${RUN_ID}"

grn_models_folder="${resources_dir}/grn_models/d0_hvgs"


reg_type=ridge
subsample=-2
num_workers=10

param_file="./params/${RUN_ID}.yaml"

grn_names=(
    'collectri'
    'negative_control'
    'positive_control'
    'pearson_corr'
    'pearson_causal'
    'portia'
    'ppcor'
    'genie3'
    'grnboost2'
    'scenic'
    'scglue'
    'celloracle'
)


# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE


append_entry() {
  cat >> $param_file << HERE
  - id: ${1}_${2}_${3}
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    reg_type: $reg_type
    method_id: ${2}-${1}
    layer: ${3}
    subsample: $subsample
    num_workers: $num_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    prediction: ${grn_models_folder}/$1.csv
    tf_all: ${resources_dir}/prior/tf_all.csv
    degree: ${2}
    noise_type: ${noise_type}
HERE
}
# Loop through grn_names and layers
layers=(scgen_pearson)
for layer in "${layers[@]}"; do
    for degree in "${degrees[@]}"; do
        for grn_name in "${grn_names[@]}"; do
            append_entry "$grn_name" "$degree" "$layer"
        done
    done
done 

# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_robustness_analysis/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

./tw launch https://github.com/openproblems-bio/run_robustness_analysis \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file ${param_file} \
  --config src/common/nextflow_helpers/labels_tw.config

