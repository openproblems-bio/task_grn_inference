workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | scenicplus.run(
      fromState: [multiomics_rna: "multiomics_rna",
              multiomics_atac: "multiomics_atac",
              temp_dir: "temp_dir",
              num_workers: "num_workers"
              ],
      toState: [prediction:"prediction", cistopic_object:"cistopic_object"]
    )

    | setState(["prediction", "cistopic_object"])

  emit:
  output_ch
}
