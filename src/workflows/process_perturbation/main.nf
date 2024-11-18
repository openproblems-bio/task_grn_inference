workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | sc_counts.run(
      fromState: [perturbation_counts: "perturbation_counts"],
      toState: [pseudobulked_data:"pseudobulked_data", pseudobulked_data_f: "pseudobulked_data_f"]
    )

    | normalization.run(
      fromState: [pseudobulked_data_f: "pseudobulked_data_f"],
      toState: [perturbation_data_n: "perturbation_data_n"]
    )

    
    | setState(["perturbation_data_n"])

  emit:
  output_ch
}
