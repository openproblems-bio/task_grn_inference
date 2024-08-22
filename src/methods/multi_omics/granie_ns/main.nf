workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | granie.run(
      fromState: [
              multiomics_rna_r: "multiomics_rna_r",
              multiomics_ata_r: "multiomics_ata_r",
              num_workers: "num_workers"
              ],
      toState: [prediction:"prediction"]
    )

    | setState(["prediction"])

  emit:
  output_ch
}
