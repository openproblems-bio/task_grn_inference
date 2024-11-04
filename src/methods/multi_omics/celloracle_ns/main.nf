workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | celloracle.run(
      fromState: [rna: "rna",
              atac: "atac",
              tf_all: "tf_all",
              temp_dir: "temp_dir",
              num_workers: "num_workers",
              max_n_links: "max_n_links"
              ],
      toState: [prediction:"prediction", base_grn: "base_grn", links: "links"]
    )

    | setState(["prediction", "base_grn", "links"])

  emit:
  output_ch
}
