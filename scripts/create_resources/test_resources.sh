#!/bin/bash
# Build the test resources (op + norman) in resources_test/task_grn_inference/
# Run after: viash ns build --parallel
set -e
cd "$(git rev-parse --show-toplevel)"

IN=resources_test/task_grn_inference/inputs
OUT=resources_test/task_grn_inference

python scripts/create_resources/subset_test_inputs.py $IN

cat > $IN/param_list.yaml <<EOF
- id: op
  dataset_id: op
  input_rna: $IN/op/rna.h5ad
  input_atac: $IN/op/atac.h5ad
  input_evaluation: $IN/op/evaluation.h5ad
- id: norman
  dataset_id: norman
  input_rna: $IN/norman/rna.h5ad
  input_evaluation: $IN/norman/evaluation.h5ad
  input_evaluation_de: $IN/norman/evaluation_de.h5ad
  input_evaluation_sc: $IN/norman/evaluation_sc.h5ad
  input_ws_background: $IN/norman/ws_background.csv
EOF

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  -c common/nextflow_helpers/labels_ci.config \
  --param_list $IN/param_list.yaml \
  --publish_dir $OUT

# shared priors (tf_all, ground truth, gene sets) are copied as-is
mkdir -p $OUT/prior
cp resources_test/grn_benchmark/prior/tf_all.csv $OUT/prior/
cp -r resources_test/grn_benchmark/prior/pathways $OUT/prior/
mkdir -p $OUT/ground_truth && cp resources_test/grn_benchmark/ground_truth/PBMC.csv $OUT/ground_truth/PBMC_unibind.csv
