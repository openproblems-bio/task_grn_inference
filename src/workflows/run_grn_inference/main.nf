workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | celloracle.run(
      fromState: [multiomics_rna: "multiomics_rna",
              multiomics_atac: "multiomics_atac",
              temp_dir: "temp_dir",
              num_workers: "num_workers"
              ],
      toState: [prediction:"prediction_celloracle"]
    )
    | setState(["prediction_celloracle"])

  emit:
  output_ch
}
