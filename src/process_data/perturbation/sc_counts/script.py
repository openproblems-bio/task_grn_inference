# !pip install sctk anndata
# !aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources_raw/ --no-sign-request

import anndata as ad 
import pandas as pd
import numpy as np
import sctk
from scipy import sparse
import scanpy as sc

par = {
    'sc_counts': 'resources/datasets_raw/perturbation_counts.h5ad',
    'perturbation_data_f': 'resources/grn-benchmark/perturbation_data.h5ad',
}

def preprocess_sc(par):
    # clean up
    sc_counts = ad.read_h5ad(par['sc_counts'])
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
def sum_by(adata: ad.AnnData, col: str) -> ad.AnnData:
    """
    Adapted from this forum post: 
    https://discourse.scverse.org/t/group-sum-rows-based-on-jobs-feature/371/4
    """
    
    assert pd.api.types.is_categorical_dtype(adata.obs[col])

    # sum `.X` entries for each unique value in `col`
    cat = adata.obs[col].values

    indicator = sparse.coo_matrix(
        (
            np.broadcast_to(True, adata.n_obs),
            (cat.codes, np.arange(adata.n_obs))
        ),
        shape=(len(cat.categories), adata.n_obs),
    )
  
    sum_adata = ad.AnnData(
        indicator @ adata.X,
        var=adata.var,
        obs=pd.DataFrame(index=cat.categories),
    )
    
    # copy over `.obs` values that have a one-to-one-mapping with `.obs[col]`
    obs_cols = adata.obs.columns
    obs_cols = list(set(adata.obs.columns) - set([col]))
    
    one_to_one_mapped_obs_cols = []
    nunique_in_col = adata.obs[col].nunique()
    for other_col in obs_cols:
        if len(adata.obs[[col, other_col]].drop_duplicates()) == nunique_in_col:
            one_to_one_mapped_obs_cols.append(other_col)

    joining_df = adata.obs[[col] + one_to_one_mapped_obs_cols].drop_duplicates().set_index(col)
    assert (sum_adata.obs.index == sum_adata.obs.join(joining_df).index).all()
    sum_adata.obs = sum_adata.obs.join(joining_df)
    sum_adata.obs.index.name = col
    sum_adata.obs = sum_adata.obs.reset_index()
    sum_adata.obs.index = sum_adata.obs.index.astype('str')

    return sum_adata
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
def filter_func(bulk_adata):
    '''Filters pseudobulked data by removing outliers compound, 
    samples with low cell counts, and genes with low coverage
    '''
    ### filter
    # samples with less than 10 cells
    bulk_adata_filtered = bulk_adata.copy()
    # toxic ones
    outliers_toxic = ['Alvocidib', 'UNII-BXU45ZH6LI', 'CGP 60474', 'BMS-387032']
    bulk_adata_filtered = bulk_adata_filtered[~bulk_adata_filtered.obs.sm_name.isin(outliers_toxic),:]
    # remove those with less than 10 cells left 
    mask_low_cell_count = bulk_adata_filtered.obs.cell_count < 10
    print(mask_low_cell_count.shape)
    bulk_adata_filtered = bulk_adata_filtered[~mask_low_cell_count]
    # remove those that have less than 2 cells types left per donor
    to_go_compounds = []
    for donor_id in bulk_adata_filtered.obs.donor_id.unique():
        adata_donor = bulk_adata_filtered[bulk_adata_filtered.obs.donor_id.eq(donor_id)]
        cell_type_n = adata_donor.obs.groupby('sm_name').size()
        to_go_compounds.append(cell_type_n[cell_type_n<=2].index.astype(str))
    to_go_compounds = np.unique(np.concatenate(to_go_compounds))
    outliers_two_celltype = ['CEP-18770 (Delanzomib)', 'IN1451', 'MLN 2238', 'Oprozomib (ONX 0912)']
    assert np.all(to_go_compounds==outliers_two_celltype)
    bulk_adata_filtered = bulk_adata_filtered[~bulk_adata_filtered.obs.sm_name.isin(to_go_compounds),:]

    # remove big class misbalance in all donors 
    outliers_misbalance_all = ['Proscillaridin A;Proscillaridin-A'] 
    bulk_adata_filtered = bulk_adata_filtered[~bulk_adata_filtered.obs.sm_name.isin(outliers_misbalance_all),:]
    # remove big class misbalance in 1 donor
    outliers_misbalance_donor_2 = ['Vorinostat']
    bulk_adata_filtered = bulk_adata_filtered[~ (bulk_adata_filtered.obs.sm_name.isin(outliers_misbalance_donor_2) & (bulk_adata_filtered.obs.donor_id=='Donor 2')),:]
    outliers_misbalance_donor_3 = ['AT13387', 'Ganetespib (STA-9090)']
    bulk_adata_filtered = bulk_adata_filtered[~ (bulk_adata_filtered.obs.sm_name.isin(outliers_misbalance_donor_3) & (bulk_adata_filtered.obs.donor_id=='Donor 3')),:]
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

    
sc_counts_f = preprocess_sc(par)
bulk_adata = pseudobulk_sum_func(sc_counts_f)
bulk_adata = pseudobulk_mean_func(bulk_adata)
bulk_adata_filtered = filter_func(bulk_adata)
bulk_adata_filtered.write(par['perturbation_data_f'])
