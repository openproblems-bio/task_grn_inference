workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | granie.run(
      fromState: [
              multiomics_rna_r: "multiomics_rna_r",
              multiomics_atac_r: "multiomics_atac_r",
              num_workers: "num_workers",
              subset: "subset"
              ],
      toState: [prediction:"prediction"]
    )

    | setState(["prediction"])

  emit:
  output_ch
}
