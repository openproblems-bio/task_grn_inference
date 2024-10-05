workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | format_data.run(
      fromState: [multiome_counts: "multiome_counts"],
      toState: [multiomics_rna: "multiomics_rna", multiomics_atac: "multiomics_atac"]
    )

    | multiome_matrix.run(
      fromState: [multiomics_rna: "multiomics_rna", multiomics_atac: "multiomics_atac"],
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
      toState: [rna_rds: "rna_rds",
                atac_rds: "atac_rds"]
    )

    | subset_hvg.run(
      fromState: [multiomics_rna: "multiomics_rna", multiomics_atac: "multiomics_atac"],
      toState: [multiomics_rna_d0_hvg: "multiomics_rna_d0_hvg", multiomics_atac_d0: "multiomics_atac_d0"]
    )

    | multiome_matrix.run(
      fromState: [multiomics_rna: "multiomics_rna_d0_hvg", multiomics_atac: "multiomics_atac_d0"],
      toState: [rna_matrix_d0: "rna_matrix",
                atac_matrix_d0: "atac_matrix",
                rna_gene_annot_d0: "rna_gene_annot",
                rna_cell_annot_d0: "rna_cell_annot",
                atac_peak_annot_d0: "atac_peak_annot",
                atac_cell_annot_d0: "atac_cell_annot"]
    )
    
    | format_resources_r.run(
      fromState: [rna_matrix: "rna_matrix_d0",
                atac_matrix: "atac_matrix_d0",
                rna_gene_annot: "rna_gene_annot_d0",
                rna_cell_annot: "rna_cell_annot_d0",
                atac_peak_annot: "atac_peak_annot_d0",
                atac_cell_annot: "atac_cell_annot_d0"],
      toState: [rna_rds_d0_hvg: "rna_rds",
                atac_rds_d0: "atac_rds"]
    )

    | setState(["multiomics_rna", "multiomics_atac", "rna_rds", "atac_rds", "multiomics_rna_d0_hvg", "multiomics_atac_d0", "rna_rds_d0_hvg", "atac_rds_d0"])

  emit:
  output_ch
}
