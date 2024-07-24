#!/bin/bash

RUN_ID="run_figr_$(date +%Y-%m-%d_%H-%M-%S)"
resources_dir="s3://openproblems-data/resources/grn"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

cat > /tmp/params.yaml << HERE
param_list:
  - id: neurips-2023-data
    de_train_h5ad: "$resources_dir/neurips-2023-data/de_train.h5ad"
    de_test_h5ad: "$resources_dir/neurips-2023-data/de_test.h5ad"
    id_map: "$resources_dir/neurips-2023-data/id_map.csv"
    layer: clipped_sign_log10_pval

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch openproblems-bio/task_perturbation_prediction \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --config src/common/nextflow_helpers/labels_tw.config
