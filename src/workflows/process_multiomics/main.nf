workflow run_wf {
  take:
  input_ch
  main:
  output_ch = input_ch
    | format_data.run(
      fromState: [multiome_counts: "multiome_counts"],
      toState: [multiomics_rna: "multiomics_rna", multiomics_atac: "multiomics_atac"]
    )
  emit:
  output_ch
}
