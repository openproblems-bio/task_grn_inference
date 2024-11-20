workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | opsca.run(
      fromState: [perturbation_counts: "perturbation_counts"],
      toState: [pseudobulked_data:"pseudobulked_data"]
    )
    
    | setState(["pseudobulked_data"])

  emit:
  output_ch
}
