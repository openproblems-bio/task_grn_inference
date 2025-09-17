dataset=$1
method="scprint"

if [ -z "$dataset" ]; then
    echo "Please provide a dataset name (e.g., 'ibd')."
    exit 1
fi

viash run src/methods/scprint/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/${dataset}_rna.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --prediction resources/results/$dataset/$dataset.$method.$method.prediction.h5ad

