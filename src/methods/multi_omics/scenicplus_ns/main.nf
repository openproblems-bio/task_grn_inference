workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | scenicplus.run(
      fromState: [
              multiomics_rna: "multiomics_rna",
              multiomics_atac: "multiomics_atac",
              temp_dir: "temp_dir",
              num_workers: "num_workers"
              
              ],
      toState: [prediction:"prediction", cell_topic:"cell_topic", scplus_mdata:"scplus_mdata"]
    )

    | setState(["prediction", "cell_topic", "scplus_mdata"])

  emit:
  output_ch
}
