workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | sc_counts.run(
      fromState: [perturbation_counts: "perturbation_counts"],
      toState: [perturbation_data_f: "perturbation_data_f"]
    )

    | normalization.run(
      fromState: [perturbation_data_f: "perturbation_data_f"],
      toState: [perturbation_data_n: "perturbation_data_n"]
    )
    
    | batch_correction_scgen.run(
      fromState: [perturbation_data_n: "perturbation_data_n"],
      toState: [perturbation_data_bc: "perturbation_data_bc"]
    )

    | batch_correction_seurat.run(
      fromState: [perturbation_data_bc: "perturbation_data_bc"],
      toState: [perturbation_data_bc: "perturbation_data_bc"]
    )

  emit:
  output_ch
}