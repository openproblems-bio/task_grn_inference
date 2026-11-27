workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | process_dataset.run(
      fromState: [
        dataset_id: "dataset_id",
        input_rna: "input_rna",
        input_atac: "input_atac",
        input_evaluation: "input_evaluation",
        input_evaluation_de: "input_evaluation_de",
        input_evaluation_sc: "input_evaluation_sc",
        input_ws_background: "input_ws_background"
      ],
      toState: [
        output_train: "output_train",
        output_train_atac: "output_train_atac",
        output_solution: "output_solution",
        output_solution_de: "output_solution_de",
        output_solution_sc: "output_solution_sc",
        output_ws_background: "output_ws_background"
      ]
    )

    // only output the files that were produced (optional inputs -> optional outputs)
    | setState { id, state ->
      def keys = ["output_train", "output_solution"]
      if (state.input_atac) keys += "output_train_atac"
      if (state.input_evaluation_de) keys += "output_solution_de"
      if (state.input_evaluation_sc) keys += "output_solution_sc"
      if (state.input_ws_background) keys += "output_ws_background"
      state.subMap(keys)
    }

  emit:
  output_ch
}
