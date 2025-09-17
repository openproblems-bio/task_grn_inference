  
workflow auto {
  findStatesTemp(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:

  // construct list of metrics
  metrics = [
    regression_1,
    regression_2,
    ws_distance,
    sem
  ]
    
  /***************************
   * RUN METRICS *
   ***************************/
  output_ch = input_ch
    | map{ id, state ->
        [id, state + ["_meta": [join_id: id]]]
      }
    

    // run all metrics
    | runEach(
      components: metrics,
      filter: { id, state, comp ->
        !state.metric_ids || state.metric_ids.contains(comp.config.name)
      },
      id: { id, state, comp ->
        id + "." + comp.config.name
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        evaluation_data: "evaluation_data",
        evaluation_data_sc: "evaluation_data_sc",
        prediction: "prediction",
        ws_distance_background: "ws_distance_background",
        subsample: "subsample",
        reg_type: "reg_type",
        apply_tf: "apply_tf",
        apply_skeleton: "apply_skeleton",
        skeleton: "skeleton",
        num_workers: "num_workers",
        regulators_consensus: "regulators_consensus",
        ws_consensus: "ws_consensus",
        tf_all: "tf_all",
        layer: "layer"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.name,
          metric_output: output.score
        ]
      }
    )

    // extract the scores
    | extract_uns_metadata.run(
      key: "extract_scores",
      fromState: [input: "metric_output"],
      toState: { id, output, state ->
        state + [
          score_uns: readYaml(output.output).uns
        ]
      }
    )

    | joinStates { ids, states ->
      assert states[0]._meta, "no _meta found in state[0]"
      // store the metric configs in a file
      def metric_configs = metrics.collect{it.config}
      def metric_configs_yaml_blob = toYamlBlob(metric_configs)
      def metric_configs_file = tempFile("metric_configs.yaml")
      metric_configs_file.write(metric_configs_yaml_blob)

      // store the scores in a file
      def score_uns = states.collect{it.score_uns}
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)

      def new_state = [
        metric_configs: metric_configs_file,
        scores: score_uns_file,
        _meta: states[0]._meta
      ]

      ["output", new_state]
    }

  emit:
  output_ch
}
