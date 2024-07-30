# download multiome from open problems and rename it
# aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/2023-09-14_kaggle_upload/2023-08-31_sc_multiome_expression_atac.h5ad ./datasets_raw --no-sign-request
# aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources/datasets_raw/ --no-sign-request
# mv ./resources/datasets_raw/2023-08-31_sc_multiome_expression_atac.h5ad ./resources/datasets_raw/multiome_counts.h5ad
# mv ./resources/datasets_raw/2023-08-31_sc_multiome_expression_atac.h5ad ./resources/datasets_raw/perturbation_counts.h5ad
# aws s3 sync ./resources/datasets_raw/ s3://openproblems-data/resources/grn/datasets_raw/ --delete

echo ">> download raw data"
aws s3 sync s3://openproblems-data/resources/grn/datasets_raw ./resources/datasets_raw/ --delete

echo ">> Create test datasets for sc perturbation and multiome" 
viash run src/process_data/test_data_counts/config.novsh.yaml

echo ">> process multiome"
viash run src/process_data/multiomics/multiome/config.vsh.yaml 

echo ">> Preprocess perturbation data"
echo "> QC, pseudobulk, and filter"
# viash run src/process_data/perturbation/sc_counts/config.vsh.yaml 
python src/process_data/perturbation/sc_counts/script.py
echo "> Normalize counts"
# viash run src/process_data/perturbation/normalization/config.vsh.yaml 
python src/process_data/perturbation/normalization/script.py
echo "> Batch correction using scgen"
# viash run src/process_data/perturbation/batch_correction_scgen/config.vsh.yaml 
python src/process_data/perturbation/batch_correction_scgen/script.py
echo "> Batch correction using seurat"
viash run src/process_data/perturbation/batch_correction_seurat/config.vsh.yaml 



# echo ">> Perturbation data: batch correction" TODO"

# echo ">> Create prior data: TODO"

# echo ">> process supp: TODO"

echo ">> Format multiomics data for R "
echo ">> Extract matrix and annotations from multiome "
viash run  src/process_data/multiomics/multiome_matrix/config.vsh.yaml 


echo ">> Construct rds files from multiomics count matrix and annotations"
viash run  src/process_data/multiomics/multiome_r/config.vsh.yaml 


echo ">> Create test resources"
echo ">> Extract matrix and annotations from multiome: test data "
viash run  src/process_data/multiomics/multiome_matrix/config.vsh.yaml --  --multiomics_rna resources_test/grn-benchmark/multiomics_rna.h5ad \
    --multiomics_atac resources_test/grn-benchmark/multiomics_atac.h5ad \
    --rna_matrix output/scRNA/X_matrix.mtx \
    --atac_matrix output/scATAC/X_matrix.mtx \
    --rna_gene_annot output/scRNA/annotation_gene.csv \
    --rna_cell_annot output/scRNA/annotation_cell.csv \
    --atac_peak_annot output/scATAC/annotation_gene.csv \
    --atac_cell_annot output/scATAC/annotation_cell.csv

echo ">> Construct rds files from multiomics count matrix and annotations: test data"
viash run  src/process_data/multiomics/multiome_r/config.vsh.yaml --  --rna_matrix output/scRNA/X_matrix.mtx \
    --atac_matrix output/scATAC/X_matrix.mtx \
    --rna_gene_annot output/scRNA/annotation_gene.csv \
    --rna_cell_annot output/scRNA/annotation_cell.csv \
    --atac_peak_annot output/scATAC/annotation_gene.csv \
    --atac_cell_annot output/scATAC/annotation_cell.csv \
    --rna_rds resources_test/grn-benchmark/multiomics_r/rna.rds \
    --atac_rds resources_test/grn-benchmark/multiomics_r/atac.rds
