#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
RUN_ID="donor_0_normalized_noncelltype_hvg"
# resources_dir="./resources/"
resources_dir="s3://openproblems-data/resources/grn"
publish_dir="${resources_dir}/results/${RUN_ID}"

reg_type=ridge
subsample=-2
num_workers=10
layer='scgen_pearson'
metric_ids="[regression_1, regression_2]"
cell_type_specific=false #for controls
normalize=true
only_hvgs=true
# method_ids="[tigress, ennet, scsgl, pidc]"
method_ids="[pearson_corr, pearson_causal, positive_control]"

param_file="./params/${RUN_ID}.yaml"

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
  - id: ${reg_type}
    metric_ids: $metric_ids
    method_ids: $method_ids
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna_0.h5ad
    multiomics_atac: ${resources_dir}/grn-benchmark/multiomics_atac_0.h5ad
    reg_type: $reg_type
    subsample: $subsample
    num_workers: $num_workers
    layer: $layer
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    tf_all: ${resources_dir}/prior/tf_all.csv
    cell_type_specific: ${cell_type_specific}
    normalize: ${normalize}
    only_hvgs: ${only_hvgs}

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_benchmark/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_benchmark/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/benchmark_donor_0_default.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config

./tw launch https://github.com/openproblems-bio/task_grn_inference \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file ${param_file} \
  --config src/common/nextflow_helpers/labels_tw.config
