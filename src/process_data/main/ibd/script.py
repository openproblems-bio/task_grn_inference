
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import os
import numpy as np
import sys
meta = {
    'resource_dir': './',
}


sys.path.append(meta['resource_dir'])
from src.process_data.helper_data import qc_perturbation, pseudobulk_sum_func, normalize_func

def add_metadata(adata):
    adata.uns['dataset_summary'] = 'IBD multiome of PBMC (scRNA+scATAC) with 38 donors, 2 disease of UC and CD, 3 perturbation of RPMI, LPS, S. enterica'
    adata.uns['dataset_description'] = 'IBD multiome of PBMC (scRNA+scATAC) with 38 donors, 2 disease of UC and CD, 3 perturbation of RPMI, LPS, S. enterica'
    adata.uns['data_reference'] = "Not published"
    adata.uns['data_url'] = ''
    adata.uns['dataset_id'] = 'ibd'
    adata.uns['dataset_name'] = 'IBD'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    return adata

def parse_donor_id(s):
    import re
    m = re.match(r'^[A-Za-z](\d+)([MF])([A-Z]+)(\d+)$', s)
    ID = int(m.group(1))
    sex = 'Male' if m.group(2)=='M' else 'Female'
    disease = m.group(3)
    condition_map = {'1':'RPMI','2':'LPS','3':'S. enterica'}
    condition = condition_map.get(m.group(4), m.group(4))
    return pd.Series([ID, sex, disease, condition])

def qc_atac(adata,
            min_fragments=1000, max_fragments=5e4,
            min_features=500, max_features=5e4,
            max_nucleosome_signal=4,
            min_tss_enrichment=2):
    """
    Filter scATAC-seq cells based on common QC metrics.
    Requires columns: nCount_ATAC, nFeature_ATAC, nucleosome_signal, TSS.enrichment
    """
    qc_mask = (
        (adata.obs['nCount_ATAC'] >= min_fragments) &
        (adata.obs['nCount_ATAC'] <= max_fragments) &
        (adata.obs['nFeature_ATAC'] >= min_features) &
        (adata.obs['nFeature_ATAC'] <= max_features) &
        (adata.obs['nucleosome_signal'] <= max_nucleosome_signal) &
        (adata.obs['TSS.enrichment'] >= min_tss_enrichment)
    )

    print(f"Keeping {qc_mask.sum()} cells out of {adata.n_obs}")
    return adata[qc_mask].copy()
map_cell_types = {
    'CD8 T cells': 'CD8T', 
    'CD4 T cells': 'CD4T', 
    'Monocytes': 'MONO', 
    'NK cells': 'NK', 
    'B cells': 'B'
}
def format_adata(adata, data_type):
    adata.obs['donorID'] = adata.obs['donorID'].astype(str)
    parsed = adata.obs['donorID'].apply(parse_donor_id)
    parsed.columns = ['ID', 'sex', 'disease', 'condition']
    adata.obs[['ID', 'sex', 'disease', 'condition']] = parsed

    adata.obs = adata.obs[['ID', 'sex', 'disease', 'condition', 'celltype', 'cellID']].rename(
        columns={'celltype':'cell_type', 'ID': 'donor_id', 'condition': 'perturbation', 'cellID': 'barcode', 'poolID': 'pool'}
        )
    stim_map = {'LPS':'LPS', 'NS':'RPMI', 'Sal':'S. enterica', 'S. salmonella': 'S. enterica'}
    adata.obs['perturbation'] = adata.obs['perturbation'].map(lambda x: stim_map.get(x, x))
    adata.obs['cell_type'] = adata.obs['cell_type'].map(map_cell_types).astype('category')
    adata.obs['condition'] = adata.obs['perturbation'].astype('str') + '_' + adata.obs['disease'].astype('str')

    adata.obs['is_control'] = adata.obs['perturbation'].apply(lambda x: x=='RPMI')
    # we need this for consistency in naming
    if data_type == 'peak':
        peak_df = adata.var.index.to_series().str.split('-', expand=True)
        peak_df.columns = ['seqname', 'start', 'end']
        # peak_df = peak_df[peak_df['seqname'].str.startswith("chr")]
        adata.var['start'] = peak_df['start'].astype(int)
        adata.var['end'] = peak_df['end'].astype(int)
        adata.var['seqname'] = peak_df['seqname'].values
        adata = adata[:, adata.var['seqname'].str.startswith("chr")].copy()
        adata.var['ranges'] = peak_df['start'].astype(str) + '-' + peak_df['end'].astype(str)
        adata.var['strand'] = '+'
        adata.var.index = adata.var['seqname'].astype(str) + ':' + adata.var['start'].astype(str) + '-' + adata.var['end'].astype(str)
    adata.obs.index.name = 'obs_id'
    adata.X = adata.X.astype(np.float32)
    return adata

data_dir = '/vol/projects/CIIM/processed/multiome/IBD/'
adata_rna = ad.read_h5ad(f'{data_dir}/rna.h5ad')
adata_atac = ad.read_h5ad(f'{data_dir}/atac.h5ad')

adata_rna = format_adata(adata_rna, data_type='rna')
adata_rna.obs['group_col'] = adata_rna.obs['perturbation'].astype(str) + '_' + \
                             adata_rna.obs['cell_type'].astype(str) + '_' + \
                             adata_rna.obs['donor_id'].astype(str) + '_' + \
                             adata_rna.obs['disease'].astype(str)
adata_rna = qc_perturbation(adata_rna, col='group_col', min_cells_per_pert=10, min_cells_per_gene=100, min_genes_per_cell=200)

adata_atac = qc_atac(adata_atac)
adata_atac = format_adata(adata_atac, data_type='peak')

adata_rna = add_metadata(adata_rna)
adata_atac = add_metadata(adata_atac)

adata_rna.obs['perturbation'] = adata_rna.obs['perturbation'].astype(str)
adata_rna.obs['cell_type'] = adata_rna.obs['cell_type'].astype(str)
adata_rna.obs['disease'] = adata_rna.obs['disease'].astype(str)
adata_rna.obs['donor_id'] = adata_rna.obs['donor_id'].astype(str)


rna_barcodes = set(adata_rna.obs['barcode'])
atac_barcodes = set(adata_atac.obs['barcode'])
common_barcodes = rna_barcodes.intersection(atac_barcodes)
print(f'Common barcodes: {len(common_barcodes)}')
adata_rna = adata_rna[adata_rna.obs['barcode'].isin(common_barcodes)].copy()
adata_atac = adata_atac[adata_atac.obs['barcode'].isin(common_barcodes)].copy()

print(f'Final RNA shape: {adata_rna.shape}')
print(f'Final ATAC shape: {adata_atac.shape}')

adata_rna = normalize_func(adata_rna)


def split_func(adata):
    group_col = 'perturbation'
    train_group = ['RPMI']
    test_groups = list(set(adata.obs[group_col].unique().tolist()) - set(train_group))

    adata_train = adata[adata.obs[group_col].isin(train_group)]
    adata_test = adata[adata.obs[group_col].isin(test_groups)]

    print(f'Shape of training set: {adata_train.shape}', f'Shape of test set: {adata_test.shape}')
    return adata_train, adata_test



adata_train_rna, adata_test_rna = split_func(adata_rna)
adata_train_atac, adata_test_atac = split_func(adata_atac)

def harmonize(adata_rna, adata_atac):
    common_barcodes = set(adata_rna.obs['barcode']).intersection(set(adata_atac.obs['barcode']))
    print(f'Common barcodes: {len(common_barcodes)}')
    adata_rna = adata_rna[adata_rna.obs['barcode'].isin(common_barcodes)].copy()
    adata_atac = adata_atac[adata_atac.obs['barcode'].isin(common_barcodes)].copy()
    return adata_rna, adata_atac
adata_train_rna, adata_train_atac= harmonize(adata_train_rna, adata_train_atac)
adata_test_rna, adata_test_atac= harmonize(adata_test_rna, adata_test_atac)

adata_test_rna_bulk = pseudobulk_sum_func(adata_test_rna, group='group_col')
adata_test_rna_bulk = normalize_func(adata_test_rna_bulk)
adata_test_rna_bulk.layers['lognorm'] = adata_test_rna_bulk.X.copy()
adata_test_rna_bulk = add_metadata(adata_test_rna_bulk)

adata_rna_bulk = pseudobulk_sum_func(adata_rna, group='group_col')
adata_rna_bulk = normalize_func(adata_rna_bulk)
adata_rna_bulk.layers['lognorm'] = adata_rna_bulk.X.copy()
adata_rna_bulk = add_metadata(adata_rna_bulk)


for disease in ['CD', 'UC']:
    adata_train_rna[adata_train_rna.obs['disease']==disease].write(f'resources/grn_benchmark/inference_data/ibd_{disease}_rna.h5ad')
    adata_train_atac[adata_train_atac.obs['disease']==disease].write(f'resources/grn_benchmark/inference_data/ibd_{disease}_atac.h5ad')
    adata_test_rna_bulk[adata_test_rna_bulk.obs['disease']==disease].write(f'resources/grn_benchmark/evaluation_data/ibd_{disease}_bulk.h5ad')
    adata_test_rna[adata_test_rna.obs['disease']==disease].write(f'resources/processed_data/ibd_{disease}_sc.h5ad')
    adata_rna_bulk[adata_rna_bulk.obs['disease']==disease].write(f'resources/extended_data/ibd_{disease}_bulk.h5ad')

print('Done')