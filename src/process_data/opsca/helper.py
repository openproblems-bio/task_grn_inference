import anndata as ad 
import pandas as pd
import numpy as np
import sctk
from scipy import sparse
import scanpy as sc

import sys


def add_metadata(adata):
    adata.uns['dataset_id'] = 'op'
    adata.uns['dataset_name'] = 'OPSCA'
    adata.uns['dataset_summary'] = 'scRNA-seq data with 146 (originally) perturbations with chemical compounds on PBMCs. Multiome data available for the control compound.'
    adata.uns['dataset_organism'] = 'human'
    adata.uns["dataset_description"] = "Novel single-cell perturbational dataset in human peripheral blood mononuclear cells (PBMCs). 144 compounds were selected from the Library of Integrated Network-Based Cellular Signatures (LINCS) Connectivity Map dataset (PMID: 29195078) and measured single-cell gene expression profiles after 24 hours of treatment. The experiment was repeated in three healthy human donors, and the compounds were selected based on diverse transcriptional signatures observed in CD34+ hematopoietic stem cells (data not released). This experiment was performed in human PBMCs because the cells are commercially available with pre-obtained consent for public release and PBMCs are a primary, disease-relevant tissue that contains multiple mature cell types (including T-cells, B-cells, myeloid cells, and NK cells) with established markers for annotation of cell types. To supplement this dataset, joint scRNA and single-cell chromatin accessibility measurements were measured from the baseline compound using the 10x Multiome assay. "
    adata.uns["data_reference"] = "@article{slazata2024benchmark,\n\ttitle = {A benchmark for prediction of transcriptomic responses to chemical perturbations across cell types},\n\tauthor = {Artur Szałata and Andrew Benz and Robrecht Cannoodt and Mauricio Cortes and Jason Fong and Sunil Kuppasani and Richard Lieberman and Tianyu Liu and Javier A. Mas-Rosario and Rico Meinl and Jalil Nourisa and Jared Tumiel and Tin M. Tunjic and Mengbo Wang and Noah Weber and Hongyu Zhao and Benedict Anchang and Fabian J Theis and Malte D Luecken and Daniel B Burkhardt},\n\tbooktitle = {The Thirty-eight Conference on Neural Information Processing Systems Datasets and Benchmarks Track},\n\tyear = {2024},\n\turl = {https://openreview.net/forum?id=WTI4RJYSVm}\n}"
    adata.uns["data_url"] = "https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP527159"
    return adata

def shorten_inference_data(par):
    # - shorten rna 
    adata_rna = ad.read_h5ad(par['op_rna'])

    n_cells = 2000
    n_peaks = 10000
    mask = adata_rna.obs.donor_id=='donor_0'
    adata_rna_s = adata_rna[mask]
    # only one chr and n_cells 
    if 'interval' not in adata_rna.var:
        raise ValueError('location is not given in rna')
    chr_mask = adata_rna_s.var.interval.str.split(':').str[0] == 'chr1'
    adata_rna_s = adata_rna_s[:n_cells, chr_mask]
   
    # - shorten atac
    adata_atac = ad.read_h5ad(par['op_atac'])
    mask = adata_atac.obs.donor_id=='donor_0'
    adata_atac_s = adata_atac[mask]
    chr_mask = adata_atac_s.var.index.str.split(':').str[0] == 'chr1'
    adata_atac_s = adata_atac_s[adata_rna_s.obs_names, chr_mask]

    total_counts = adata_atac_s.X.sum(axis=0).A1  # .A1 converts the sparse matrix to a dense array
    peaks_df = pd.DataFrame({
        'peak': adata_atac_s.var.index,
        'total_counts': total_counts
    })
    top_peaks = peaks_df.nlargest(n_peaks, 'total_counts')['peak']
    adata_atac_s = adata_atac_s[:, top_peaks]
    adata_atac_s.var = adata_atac_s.var.loc[top_peaks]

    print(adata_rna_s)
    print(adata_atac_s)
    assert adata_rna_s.shape[0] > 0, 'no cells in rna'
    assert adata_atac_s.shape[0] > 0, 'no cells in atac'
    assert adata_rna_s.shape[1] > 0, 'no genes in rna'
    
    adata_rna_s.write(par['op_rna_test'])
    adata_atac_s.write(par['op_atac_test'])
def shorten_evaluation_data(par):
    evaluation_data = ad.read_h5ad(par['op_perturbation_bulk'])
    evaluation_data = evaluation_data[evaluation_data.obs.sample(2000).index, evaluation_data.var.sample(2000).index]
    evaluation_data.write(par['op_perturbation_bulk_test'])

def preprocess_sc(par):
    # clean up
    sc_counts = ad.read_h5ad(par['op_perturbation_raw'])
    sc_counts.obs = sc_counts.obs[['well', 'row', 'col', 'plate_name', 'cell_type', 'donor_id', 'sm_name']]
    sc_counts.X = sc_counts.layers['counts']
    del sc_counts.layers 
    del sc_counts.obsm 
    sc_counts.var_names_make_unique()
    # merge cell types
    CELL_TYPES = ['NK cells', 'T cells CD4+', 'T cells CD8+', 'T regulatory cells', 'B cells', 'Myeloid cells']
    T_cell_types = ['T regulatory cells', 'T cells CD8+', 'T cells CD4+']
    cell_type_map = {cell_type: 'T cells' if cell_type in T_cell_types else cell_type for cell_type in CELL_TYPES}
    sc_counts.obs['cell_type'] = sc_counts.obs['cell_type'].map(cell_type_map)
    sc_counts.obs['cell_type'].unique()

    # qc 
    sctk.calculate_qc(sc_counts)
    sctk.cellwise_qc(sc_counts)

    # filtering
    # cell wise
    filter_percent_hb = sc_counts.obs.percent_hb>.2
    filter_percent_hb.sum()
    # gene wise
    plates = sc_counts.obs['plate_name'].unique()

    # Step 2: Initialize a DataFrame to store counts
    gene_counts_per_plate = pd.DataFrame(index=sc_counts.var_names, columns=plates, dtype=int)

    # Step 3: Iterate over each plate and calculate expression counts
    for plate in plates:
        # Subset the AnnData object for the current plate
        subset = sc_counts[sc_counts.obs['plate_name'] == plate]

        # Calculate expression counts (genes x cells > 0)
        expressed_genes = (subset.X > 0).sum(axis=0)

        # Check if the result needs conversion from sparse matrix format
        if isinstance(expressed_genes, np.matrix):
            expressed_genes = np.array(expressed_genes).flatten()

        # Store the counts in the DataFrame
        gene_counts_per_plate[plate] = expressed_genes

    # Step 4: Aggregate counts across plates (max or sum based on the requirement)
    # We use `max` here to find if any gene meets the criteria in at least one plate
    max_counts = gene_counts_per_plate.max(axis=1)

    # Step 5: Create a mask for genes to keep (genes expressed in at least 100 cells in any plate)
    genes_to_keep = max_counts >= 100
    print('retained genes:', genes_to_keep.sum())
    # actual filtering
    sc_counts = sc_counts[(~filter_percent_hb), genes_to_keep]
    # clean
    sc_counts.obs = sc_counts.obs[['cell_type', 'sm_name', 'donor_id', 'row', 'plate_name', 'well']]
    sc_counts.var = sc_counts.var[[]]

    del sc_counts.obsm
    del sc_counts.uns
    return sc_counts


def pseudobulk_mean_func(bulk_adata):
    bulk_adata.layers['counts'] = bulk_adata.X.copy()
    rows_adj = []
    for i, row in enumerate(bulk_adata.X):
        count = bulk_adata.obs.cell_count[i]
        rows_adj.append(row/count)

    bulk_adata.layers['n_counts'] = np.asarray(rows_adj)

    return bulk_adata
def filter_func(bulk_adata, cell_counts_t):
    '''Filters pseudobulked data by removing outliers compound, 
    samples with low cell counts, and genes with low coverage
    '''
    ### filter
    # samples with less than 10 cells
    bulk_adata_filtered = bulk_adata.copy()
    bulk_adata_filtered.obs['donor_id'] = bulk_adata_filtered.obs.donor_id.map({'Donor 1': 'donor_0', 'Donor 2': 'donor_1', 'Donor 3': 'donor_2'})
    # toxic ones
    outliers_toxic = ['Alvocidib', 'UNII-BXU45ZH6LI', 'CGP 60474', 'BMS-387032']
    bulk_adata_filtered = bulk_adata_filtered[~bulk_adata_filtered.obs.sm_name.isin(outliers_toxic),:]
    # remove those with less than 10 cells left 

    mask_low_cell_count = bulk_adata_filtered.obs.cell_count < cell_counts_t
    bulk_adata_filtered = bulk_adata_filtered[~mask_low_cell_count]
    
    # remove those that have less than 2 cells types left per donor
    to_go_compounds = []
    for donor_id in bulk_adata_filtered.obs.donor_id.unique():
        adata_donor = bulk_adata_filtered[bulk_adata_filtered.obs.donor_id.eq(donor_id)]
        cell_type_n = adata_donor.obs.groupby('sm_name').size()
        to_go_compounds.append(cell_type_n[cell_type_n<=2].index.astype(str))
    to_go_compounds = np.unique(np.concatenate(to_go_compounds))
    outliers_two_celltype = ['CEP-18770 (Delanzomib)', 'IN1451', 'MLN 2238', 'Oprozomib (ONX 0912)']
    # assert np.all(to_go_compounds==outliers_two_celltype)
    bulk_adata_filtered = bulk_adata_filtered[~bulk_adata_filtered.obs.sm_name.isin(to_go_compounds),:]

    # remove big class misbalance in all donors 
    outliers_misbalance_all = ['Proscillaridin A;Proscillaridin-A'] 
    bulk_adata_filtered = bulk_adata_filtered[~bulk_adata_filtered.obs.sm_name.isin(outliers_misbalance_all),:]
    # remove big class misbalance in 1 donor
    outliers_misbalance_donor_2 = ['Vorinostat']
    bulk_adata_filtered = bulk_adata_filtered[~ (bulk_adata_filtered.obs.sm_name.isin(outliers_misbalance_donor_2) & (bulk_adata_filtered.obs.donor_id=='donor_1')),:]
    outliers_misbalance_donor_3 = ['AT13387', 'Ganetespib (STA-9090)']
    bulk_adata_filtered = bulk_adata_filtered[~ (bulk_adata_filtered.obs.sm_name.isin(outliers_misbalance_donor_3) & (bulk_adata_filtered.obs.donor_id=='donor_2')),:]
    print(f"number of initial samples: {len(bulk_adata)}, number of samples after filtering: {len(bulk_adata_filtered)}")
    # low gene coverage
    mask_to_go_genes = ((bulk_adata_filtered.X == 0).sum(axis=0)/bulk_adata_filtered.shape[0])>0.7
    print('number of removed genes:', mask_to_go_genes.sum())
    bulk_adata_filtered = bulk_adata_filtered[:,~mask_to_go_genes] 
    bulk_adata_filtered.obs.drop(columns=['plate_well_cell_type'], inplace=True)
    # for the sake of seurat
    for key in ['cell_type','plate_name']:
        bulk_adata_filtered.obs[key] = bulk_adata_filtered.obs[key].astype(str)
    return bulk_adata_filtered
def normalize_func(adata):
    # sc.pp.normalize_total(bulk_adata_c)
    # sc.pp.log1p(bulk_adata_c)

    adata.layers['X_norm'] = sc.experimental.pp.normalize_pearson_residuals(adata, layer='counts', inplace=False)['X']
    X_norm = sc.pp.normalize_total(adata, layer='counts',inplace=False)['X']
    X_norm = sc.pp.log1p(X_norm, copy=False)
    adata.layers['lognorm'] = X_norm
    
    return adata