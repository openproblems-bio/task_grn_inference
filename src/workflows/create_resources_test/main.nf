workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | create_test_data.run(
      fromState: [multiomics_rna: "multiomics_rna", 
                  multiomics_atac: "multiomics_atac",
                  perturbation_data: "perturbation_data"],
      toState: [multiomics_rna_test: "multiomics_rna_test", 
                multiomics_atac_test: "multiomics_atac_test",
                perturbation_data_test: "perturbation_data_test"]
    )


    | multiome_matrix.run(
      fromState: [multiomics_rna: "multiomics_rna_test", multiomics_atac: "multiomics_atac_test"],
      toState: [rna_matrix: "rna_matrix",
                atac_matrix: "atac_matrix",
                rna_gene_annot: "rna_gene_annot",
                rna_cell_annot: "rna_cell_annot",
                atac_peak_annot: "atac_peak_annot",
                atac_cell_annot: "atac_cell_annot"]
    )
    
    | format_resources_r.run(
      fromState: [rna_matrix: "rna_matrix",
                atac_matrix: "atac_matrix",
                rna_gene_annot: "rna_gene_annot",
                rna_cell_annot: "rna_cell_annot",
                atac_peak_annot: "atac_peak_annot",
                atac_cell_annot: "atac_cell_annot"],
      toState: [rna_rds: "rna_rds_test",
                atac_rds: "atac_rds_test"]
    )

    | setState(["perturbation_data_test", "multiomics_rna_test", "multiomics_atac_test", "rna_rds_test", "atac_rds_test"])

  emit:
  output_ch
}
