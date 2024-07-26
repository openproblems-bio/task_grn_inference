# download multiome and rename it
echo ">> download raw data"
# aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/2023-09-14_kaggle_upload/2023-08-31_sc_multiome_expression_atac.h5ad ./datasets_raw --no-sign-request
# aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./datasets_raw/ --no-sign-request
# mv ./datasets_raw/2023-08-31_sc_multiome_expression_atac.h5ad ./datasets_raw/multiome.h5ad
# aws s3 sync ./datasets_raw/ s3://openproblems-data/resources/datasets_raw/grn 

aws s3 sync s3://openproblems-data/resources/datasets_raw/grn ./datasets_raw/

echo ">> process multiome"
# viash run src/process_data/multiome/config.vsh.yaml -- --multiome_counts datasets_raw/multiome.h5ad --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad --multiomics_atac resources/grn-benchmark/multiomics_atac.h5ad
# python src/process_data/multiome/script.py
# echo ">> process perturbation data"
# viash run src/process_data/sc_counts/config.vsh.yaml -- --sc_counts datasets_raw/sc_counts.h5ad --perturbation_data resources/grn-benchmark/perturbation_data.h5ad 

# echo ">> Perturbation data: batch correction" TODO"

# echo ">> Create prior data: TODO"

# echo ">> process supp: TODO"

# echo ">> Extract matrix and annotations from multiome "
# viash run  src/process_data/multiome_matrix/config.novsh.yaml --  --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
#     --multiomics_atac resources/grn-benchmark/multiomics_atac.h5ad \
#     --rna_matrix output/scRNA/X_matrix.mtx \
#     --atac_matrix output/scATAC/X_matrix.mtx \
#     --rna_gene_annot output/scRNA/annotation_gene.csv \
#     --rna_cell_annot output/scRNA/annotation_cell.csv \
#     --atac_peak_annot output/scATAC/annotation_gene.csv \
#     --atac_cell_annot output/scATAC/annotation_cell.csv

# python src/process_data/multiome_matrix/script.py --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad \
#     --multiomics_atac resources/grn-benchmark/multiomics_atac.h5ad \
#     --rna_matrix output/scRNA/X_matrix.mtx \
#     --atac_matrix output/scATAC/X_matrix.mtx \
#     --rna_gene_annot output/scRNA/annotation_gene.csv \
#     --rna_cell_annot output/scRNA/annotation_cell.csv \
#     --atac_peak_annot output/scATAC/annotation_gene.csv \
#     --atac_cell_annot output/scATAC/annotation_cell.csv

# echo ">> Construct rds files from multiomics count matrix and annotations"
# Rscript src/process_data/multiome_r/script.R


echo ">> create test resources "
python src/process_data/test_data/script.py