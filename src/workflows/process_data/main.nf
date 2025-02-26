workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | op_multiome.run(
      fromState: [op_multiome: "op_multiome"],
      toState: [op_rna: "op_rna", op_atac: "op_atac"]
    )

    | op_perturbation.run(
      fromState: [op_perturbation_raw: "op_perturbation_raw"],
      toState: [op_perturbation_bulk: "op_perturbation_bulk"]
    )



    | setState(["op_rna", "op_atac", "op_perturbation_bulk"])

  emit:
  output_ch
}
