
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
    regression_1
  ]
    
  /***************************
   * RUN METRICS *
   ***************************/
  score_ch = input_ch
    | map{ id, state ->
        [id, state + ["_meta": [join_id: id]]]
      }
    // run all metrics
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "." + comp.config.functionality.name
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        layer: "layer", 
        prediction: "prediction",
        subsample: "subsample",
        reg_type: "reg_type"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.functionality.name,
          metric_output: output.score
        ]
      }
    )

  output_ch = score_ch

    // extract the scores
    | extract_metadata.run(
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

      def task_info_file = meta.resources_dir.resolve("task_info.yaml")

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

    // merge all of the output data 
    | joinStates{ ids, states ->
      def mergedStates = states.inject([:]) { acc, m -> acc + m }
      [ids[0], mergedStates]
    }

  emit:
  output_ch
}
