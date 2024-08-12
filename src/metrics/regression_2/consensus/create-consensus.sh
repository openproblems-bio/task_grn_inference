viash run src/metrics/regression_2/consensus/config.novsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
     --output resources/grn-benchmark/consensus-num-regulators.json \
     --grn_folder resources/grn_models/ \
     --grns ananse.csv,celloracle.csv,figr.csv,granie.csv,scenicplus.csv,scglue.csv