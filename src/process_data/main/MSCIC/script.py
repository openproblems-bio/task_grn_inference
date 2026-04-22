"""
MSCIC (Multimodal Single-Cell Integration Challenge) data processing script.

Source: GSE194122 — NeurIPS 2021 Open Problems Multimodal Single-Cell Integration Challenge
Tissue: Bone Marrow Mononuclear Cells (BMMC), healthy
Protocol: 10x Multiome (RNA + ATAC)
Donors: 10 unique donors across 4 sites

Split strategy:
  - Evaluation: donor 15078 @ site2 + site4, donor 18303 @ site3
                (multi-site replicates → drives replicate_consistency metric)
  - Inference:  all remaining cells (~52k cells, 8 donor-site groups)

Outputs:
  - resources/grn_benchmark/inference_data/MSCIC_rna.h5ad
  - resources/grn_benchmark/inference_data/MSCIC_atac.h5ad
  - resources/grn_benchmark/evaluation_data/MSCIC_rna_all.h5ad  (pseudobulked eval RNA)
"""

import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

sys.path.insert(0, 'src/process_data')
from helper_data import normalize_func

# ── paths ──────────────────────────────────────────────────────────────────────
par = {
    'raw': 'resources/datasets_raw/MSCIC.h5ad',
    'MSCIC_rna':       'resources/grn_benchmark/inference_data/MSCIC_rna.h5ad',
    'MSCIC_atac':      'resources/grn_benchmark/inference_data/MSCIC_atac.h5ad',
    'MSCIC_sc':        'resources/grn_benchmark/evaluation_data/MSCIC_sc.h5ad',
    'MSCIC_rna_all':      'resources/extended_data/MSCIC_rna_all.h5ad',
}

# ── evaluation mask ────────────────────────────────────────────────────────────
EVAL_DONOR_SITE = {
    '15078': {'site2', 'site4'},
    '18303': {'site3'},
}

# ── cell type harmonisation ────────────────────────────────────────────────────
CELL_TYPE_MAP = {
    'CD8+ T':           'CD8T',
    'CD8+ T naive':     'CD8T',
    'CD4+ T naive':     'CD4T',
    'CD4+ T activated': 'CD4T',
    'CD14+ Mono':       'Mono',
    'CD16+ Mono':       'Mono',
    'NK':               'NK',
    'Naive CD20+ B':    'B',
    'Transitional B':   'B',
    'B1 B':             'B',
    'Plasma cell':      'B',
    'Erythroblast':     'Erythroid',
    'Proerythroblast':  'Erythroid',
    'Normoblast':       'Erythroid',
    'HSC':              'HSC',
    'Lymph prog':       'Progenitor',
    'G/M prog':         'Progenitor',
    'MK/E prog':        'Progenitor',
    'ID2-hi myeloid prog': 'Progenitor',
    'pDC':              'DC',
    'cDC2':             'DC',
    'ILC':              'ILC',
}


def add_metadata(adata, split_label):
    adata.uns['dataset_id']      = 'MSCIC'
    adata.uns['dataset_name']    = 'MSCIC'
    adata.uns['dataset_summary'] = (
        'NeurIPS 2021 Multimodal Single-Cell Integration Challenge. '
        '10x Multiome (RNA+ATAC) of healthy bone marrow mononuclear cells (BMMCs) '
        'from 10 donors across 4 sites.'
    )
    adata.uns['dataset_description'] = (
        'GSE194122: healthy BMMC single-cell multiome dataset generated for the '
        'NeurIPS 2021 Open Problems in Single-Cell Analysis competition. '
        'Contains paired gene expression and chromatin accessibility for ~69k cells '
        'from 10 unique donors measured at up to 3 independent sites.'
    )
    adata.uns['dataset_organism'] = 'human'
    adata.uns['data_reference']   = (
        '@article{luecken2021sandbox,\n'
        '\ttitle={A sandbox for prediction and integration of DNA, RNA, and proteins '
        'in single cells},\n'
        '\tauthor={Luecken, Malte D and others},\n'
        '\tbooktitle={Thirty-fifth Conference on Neural Information Processing Systems '
        'Datasets and Benchmarks Track},\n'
        '\tyear={2021}\n}'
    )
    adata.uns['data_url']        = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122'
    adata.uns['normalization_id'] = 'lognorm'
    adata.uns['split']            = split_label
    return adata


def build_eval_mask(obs):
    """Return boolean mask for evaluation cells."""
    mask = pd.Series(False, index=obs.index)
    for donor, sites in EVAL_DONOR_SITE.items():
        mask |= (obs['DonorID'].astype(str) == donor) & (obs['Site'].isin(sites))
    return mask


def harmonise_obs(obs, donor_map):
    """Keep and rename only the columns we need."""
    out = pd.DataFrame(index=obs.index)
    out['donor_id']  = obs['DonorID'].astype(str).map(donor_map)
    out['cell_type_orig'] = obs['cell_type'].map(CELL_TYPE_MAP).fillna(obs['cell_type'])
    out['cell_type'] = 'pbmc'
    out['site']      = obs['Site'].astype(str)
    out['donor_age'] = obs['DonorAge'].values
    out['donor_gender'] = obs['DonorGender'].astype(str)
    # replicate_consistency grouping column (no perturbation → use cell_type)
    out['perturbation'] = 'healthy'
    return out


def process_rna(rna_raw, donor_map):
    """Format RNA AnnData: raw counts in X, lognorm in layers['lognorm']."""
    rna = rna_raw.copy()

    # ensure X holds raw integer counts
    X = rna.X
    if issparse(X):
        X = X.toarray()
    rna.X = csr_matrix(X.astype(np.float32))
    rna.layers['counts'] = rna.X.copy()

    # keep minimal var info
    rna.var = rna.var[['gene_id']].copy() if 'gene_id' in rna.var.columns else rna.var[[]]
    rna.var.index.name = 'gene'

    # clean obs
    rna.obs = harmonise_obs(rna_raw.obs, donor_map)
    rna.obs.index.name = 'obs_id'

    # normalise: adds layers['lognorm'], resets X to raw counts
    rna = normalize_func(rna, log_norm=True, pearson_residual=False)

    # drop obsm/obsp/varm to keep file small
    rna.obsm  = {}
    rna.obsp  = {}
    rna.varm  = {}
    rna.uns   = {}

    return rna


def process_atac(atac_raw, donor_map):
    """Format ATAC AnnData: raw counts in X, peak coordinates in var."""
    atac = atac_raw.copy()

    X = atac.X
    if issparse(X):
        X = X.toarray()
    atac.X = csr_matrix(X.astype(np.float32))

    # parse peak coordinates from var index (chrN-start-end)
    # rename to chrN:start-end format (required by scenicplus and other ATAC methods)
    atac.var = atac.var[[]]
    peak_series = pd.Series(atac.var.index, index=atac.var.index)
    peak_parts  = peak_series.str.split('-', n=2, expand=True)
    atac.var['seqname'] = peak_parts.iloc[:, 0].values
    atac.var['ranges']  = (peak_parts.iloc[:, 1] + '-' + peak_parts.iloc[:, 2]).values
    atac.var['strand']  = '+'
    # rename index to chr:start-end format for scenicplus compatibility
    new_index = (peak_parts.iloc[:, 0] + ':' +
                 peak_parts.iloc[:, 1] + '-' +
                 peak_parts.iloc[:, 2])
    atac.var.index = new_index.values
    atac.var.index.name = 'location'
    atac = atac[:, atac.var['seqname'].str.startswith('chr')]

    atac.obs = harmonise_obs(atac_raw.obs, donor_map)
    atac.obs.index.name = 'obs_id'

    atac.obsm = {}
    atac.obsp = {}
    atac.varm = {}
    atac.uns  = {}

    return atac


def pseudobulk_rna(rna, groupby_cols):
    """Sum raw counts per (donor × cell_type), then lognorm."""
    rna.obs['_group'] = rna.obs[groupby_cols].astype(str).agg('_'.join, axis=1)
    groups = rna.obs['_group'].unique()

    X_bulk = []
    obs_rows = []
    for g in groups:
        mask = rna.obs['_group'] == g
        sub = rna[mask]
        x = sub.X
        if issparse(x):
            x = x.toarray()
        X_bulk.append(x.sum(axis=0))
        row = sub.obs.iloc[0][groupby_cols].to_dict()
        row['cell_count'] = mask.sum()
        obs_rows.append(row)

    bulk = ad.AnnData(
        X=csr_matrix(np.vstack(X_bulk).astype(np.float32)),
        obs=pd.DataFrame(obs_rows),
        var=rna.var.copy(),
    )
    bulk.obs.index = [f'sample_{i}' for i in range(len(bulk.obs))]
    bulk.obs.index.name = 'obs_id'
    bulk.var.index.name = 'gene'
    bulk.layers['counts'] = bulk.X.copy()

    # lognorm
    sc.pp.normalize_total(bulk, target_sum=1e4, layer='counts')
    sc.pp.log1p(bulk)
    bulk.layers['lognorm'] = csr_matrix(bulk.X.copy())
    bulk.X = bulk.layers['counts'].copy()
    del bulk.layers['counts']

    return bulk


def main():
    print('Loading MSCIC raw data...', flush=True)
    adata = ad.read_h5ad(par['raw'], backed='r')
    print(f'  Shape: {adata.shape}', flush=True)

    # ── split features into RNA and ATAC ──────────────────────────────────────
    print('Splitting modalities...', flush=True)
    rna_mask_var  = adata.var['feature_types'] == 'GEX'
    atac_mask_var = adata.var['feature_types'] == 'ATAC'

    # load into memory for subsetting
    adata = adata.to_memory()

    rna_all  = adata[:, rna_mask_var].copy()
    atac_all = adata[:, atac_mask_var].copy()
    print(f'  RNA:  {rna_all.shape}', flush=True)
    print(f'  ATAC: {atac_all.shape}', flush=True)

    # ── build evaluation mask ─────────────────────────────────────────────────
    eval_mask = build_eval_mask(adata.obs)
    print(f'  Evaluation cells: {eval_mask.sum()}, Inference cells: {(~eval_mask).sum()}', flush=True)

    # ── anonymise donor IDs ───────────────────────────────────────────────────
    all_donors   = sorted(adata.obs['DonorID'].astype(str).unique())
    donor_map    = {d: f'donor_{i}' for i, d in enumerate(all_donors)}
    print(f'  Donor map: {donor_map}', flush=True)

    # ── split ─────────────────────────────────────────────────────────────────
    rna_infer  = rna_all[~eval_mask].copy()
    rna_eval   = rna_all[eval_mask].copy()
    atac_infer = atac_all[~eval_mask].copy()

    # ── downsample inference to ~25k cells ────────────────────────────────────
    # Always keep all cells from eval donors (multi-site donors needed for
    # replicate_consistency). Add remaining donors smallest-first until >25k.
    infer_obs = rna_infer.obs
    eval_donor_ids = set(EVAL_DONOR_SITE.keys())
    required_mask = infer_obs['DonorID'].astype(str).isin(eval_donor_ids)
    required_idx  = infer_obs.index[required_mask].tolist()
    other_idx_by_donor = (
        infer_obs[~required_mask]
        .groupby('DonorID', observed=True)
        .apply(lambda g: g.index.tolist())
        .sort_values(key=lambda s: s.map(len))  # smallest donor first
    )
    selected_idx = required_idx.copy()
    for donor_cells in other_idx_by_donor:
        selected_idx.extend(donor_cells)
        if len(selected_idx) >= 25000:
            break
    print(f'  Downsampled inference: {len(selected_idx)} cells '
          f'({len(other_idx_by_donor) + 1 - (len(selected_idx) < 25000)} donors)', flush=True)
    rna_infer  = rna_infer[selected_idx].copy()
    atac_infer = atac_infer[selected_idx].copy()

    # ── process and format ────────────────────────────────────────────────────
    print('Processing inference RNA...', flush=True)
    rna_infer = process_rna(rna_infer, donor_map)
    rna_infer = add_metadata(rna_infer, 'inference')
    print(rna_infer)

    print('Processing inference ATAC...', flush=True)
    atac_infer = process_atac(atac_infer, donor_map)
    atac_infer = add_metadata(atac_infer, 'inference')
    print(atac_infer)

    print('Processing evaluation RNA (single-cell)...', flush=True)
    rna_eval = process_rna(rna_eval, donor_map)
    rna_eval = add_metadata(rna_eval, 'evaluation')
    print(rna_eval)

    # full RNA (inference + eval) for positive_control
    print('Processing full RNA (all cells)...', flush=True)
    rna_full = process_rna(rna_all, donor_map)
    rna_full = add_metadata(rna_full, 'full')
    print(rna_full)

    # ── write ─────────────────────────────────────────────────────────────────
    print(f'Writing {par["MSCIC_rna"]}...', flush=True)
    rna_infer.write_h5ad(par['MSCIC_rna'])

    print(f'Writing {par["MSCIC_atac"]}...', flush=True)
    atac_infer.write_h5ad(par['MSCIC_atac'])

    print(f'Writing {par["MSCIC_sc"]}...', flush=True)
    rna_eval.write_h5ad(par['MSCIC_sc'])

    print(f'Writing {par["MSCIC_rna_all"]}...', flush=True)
    rna_full.write_h5ad(par['MSCIC_rna_all'])

    print('\nDone.', flush=True)
    print(f'  Inference RNA:  {rna_infer.shape}')
    print(f'  Inference ATAC: {atac_infer.shape}')
    print(f'  Evaluation SC:  {rna_eval.shape}')
    print(f'  Full RNA:       {rna_full.shape}')


if __name__ == '__main__':
    main()
