# !pip install sctk anndata
# !aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources_raw/ --no-sign-request

import anndata as ad 
import pandas as pd
import numpy as np
import sctk
from scipy import sparse
import scanpy as sc
from scipy.sparse import csr_matrix
import sys

## VIASH START
par = {
    'perturbation_counts': 'resources/datasets_raw/op_perturbation_sc_counts.h5ad',
    'perturbation_bulk': 'resources/extended_data/op_perturbation_bulk.h5ad',
    'evaluation_data': 'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
}
## VIASH END
meta = { 
    'resources': 'src/utils/'
}   

sys.path.append(meta['resources'])

from util import sum_by

def preprocess_sc(par):
    # clean up
    sc_counts = ad.read_h5ad(par['perturbation_counts'])
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

def pseudobulk_sum_func(sc_counts):
    # pseudobulk
    #group cell types per well
    sc_counts.obs['plate_well_cell_type'] = sc_counts.obs['plate_name'].astype('str') \
        + '_' + sc_counts.obs['well'].astype('str') \
        + '_' + sc_counts.obs['cell_type'].astype('str')
    sc_counts.obs['plate_well_cell_type'] = sc_counts.obs['plate_well_cell_type'].astype('category')
    bulk_adata = sum_by(sc_counts, 'plate_well_cell_type')
    bulk_adata.obs['cell_count'] = sc_counts.obs.groupby('plate_well_cell_type').size().values
    bulk_adata.X = np.array(bulk_adata.X.todense())

    print('ratio of missingness' , (bulk_adata.X==0).sum()/bulk_adata.X.size)
    bulk_adata.var = bulk_adata.var.reset_index()
    bulk_adata.var.set_index('index', inplace=True)

    bulk_adata.X = np.nan_to_num(bulk_adata.X, nan=0)
    return bulk_adata
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

if __name__ == '__main__':
    
    cell_counts_t = 10
        
    sc_counts_f = preprocess_sc(par)
    bulk_adata = pseudobulk_sum_func(sc_counts_f)
    bulk_adata = pseudobulk_mean_func(bulk_adata)
    bulk_adata = filter_func(bulk_adata, cell_counts_t)


    bulk_adata.obs = bulk_adata.obs.rename(columns={'sm_name':'perturbation'})

    bulk_adata = normalize_func(bulk_adata)

    del bulk_adata.layers['n_counts']

    bulk_adata.X = csr_matrix(bulk_adata.X)
    bulk_adata.layers['counts'] = csr_matrix(bulk_adata.layers['counts'])

    bulk_adata.obs['is_control'] = bulk_adata.obs['perturbation'].isin(['Dimethyl Sulfoxide'])
    bulk_adata.obs['is_positive_control'] = bulk_adata.obs['perturbation'].isin(['Dabrafenib', 'Belinostat'])

    bulk_adata.uns['dataset_id'] = 'op'
    bulk_adata.write(par['perturbation_bulk'])
    bulk_adata.write(par['evaluation_data'])