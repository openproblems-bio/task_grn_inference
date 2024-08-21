workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | granie.run(
      fromState: [
              multiomics_rna: "multiomics_rna",
              multiomics_atac: "multiomics_atac",
              num_workers: "num_workers"
              ],
      toState: [prediction:"prediction"]
    )

    | setState(["prediction"])

  emit:
  output_ch
}
