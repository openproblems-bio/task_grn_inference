"""
SoundLife sc CD4T-cell processing — Scenario 1a.

Temporal split: same donors, same condition (day 0, pre-vaccination baseline),
across two consecutive flu years.

  Inference:  CD4T cells, day=0, year=1  (10 donors)
  Evaluation: CD4T cells, day=0, year=2  (same 10 donors)

Goal: evaluate whether a GRN inferred from resting CD4T cells in year 1
generalises to the same donors' resting state one year later.

Source: /vol/projects/jnourisa/hiara/datasets_old/sc/soundlife.h5ad
Outputs:
  resources/grn_benchmark/inference_data/soundlife_rna_s1.h5ad
  resources/grn_benchmark/evaluation_data/soundlife_sc_s1.h5ad
"""

import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix, issparse

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src', 'process_data'))
from helper_data import normalize_func

RAW   = '/vol/projects/jnourisa/hiara/datasets_old/sc/soundlife.h5ad'
N_DONORS = 10
CELLS_PER_DONOR = 2500

par = {
    'inference':  os.path.join(ROOT, 'resources/grn_benchmark/inference_data/soundlife_rna.h5ad'),
    'evaluation': os.path.join(ROOT, 'resources/grn_benchmark/evaluation_data/soundlife_sc.h5ad'),
    'rna_all':    os.path.join(ROOT, 'resources/extended_data/soundlife_rna_all.h5ad'),
}


def downsample_per_donor(adata, cells_per_donor=CELLS_PER_DONOR, seed=0):
    rng = np.random.default_rng(seed)
    keep = []
    for donor in adata.obs['donor_id'].unique():
        idx = adata.obs.index[adata.obs['donor_id'] == donor].tolist()
        n = min(cells_per_donor, len(idx))
        keep.extend(rng.choice(idx, size=n, replace=False).tolist())
    return adata[keep].copy()


def select_donors(obs, infer_mask, eval_mask, n=N_DONORS, seed=0):
    """Return n donor_ids present in both inference and evaluation."""
    infer_donors = set(obs.loc[infer_mask, 'donor_id'].unique())
    eval_donors  = set(obs.loc[eval_mask,  'donor_id'].unique())
    shared = sorted(infer_donors & eval_donors)
    rng = np.random.default_rng(seed)
    return list(rng.choice(shared, size=min(n, len(shared)), replace=False))


def harmonise_obs(obs, donor_map):
    out = pd.DataFrame(index=obs.index)
    out['donor_id']      = obs['donor_id'].map(donor_map)
    out['cell_type']     = 'CD4T'
    out['cell_type_orig'] = obs['Sub_CT'].astype(str)
    out['perturbation']  = 'pre_vaccination'
    out['day']           = obs['day'].astype(str)
    out['year']          = obs['year'].astype(str)
    out['visitName']     = obs['visitName'].astype(str)
    out['donor_age']     = obs['sample.subjectAgeAtDraw'].astype(str)
    out['donor_sex']     = obs['subject.biologicalSex'].astype(str)
    out.index.name = 'obs_id'
    return out


def add_metadata(adata, split_label):
    adata.uns['dataset_id']       = 'soundlife'
    adata.uns['dataset_name']     = 'SoundLife'
    adata.uns['dataset_summary']  = (
        'SoundLife longitudinal flu-vaccination cohort. '
        'CD4T cells from healthy adults sampled at day 0 (pre-vaccination baseline) '
        'across two flu seasons.'
    )
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    adata.uns['split']            = split_label
    adata.uns['scenario']         = '1a'
    adata.uns['scenario_description'] = (
        'Inference: day=0, year=1. '
        'Evaluation: day=0, year=2 (same donors, one year later).'
    )
    return adata


def process(adata_subset, donor_map):
    adata = adata_subset.copy()
    X = adata.X
    if issparse(X):
        X = X.toarray()
    adata.X = csr_matrix(X.astype(np.float32))
    adata.obs = harmonise_obs(adata_subset.obs, donor_map)
    adata.var = adata.var[[]] if len(adata.var.columns) else adata.var
    adata.var.index.name = 'gene'
    adata.obsm = {}
    adata.obsp = {}
    adata.varm = {}
    adata.uns  = {}
    adata = normalize_func(adata, log_norm=True)
    return adata


def main():
    print('Loading backed...', flush=True)
    raw = sc.read_h5ad(RAW, backed='r')
    print(f'  {raw.shape}', flush=True)

    obs = raw.obs.copy()

    # ── define splits ─────────────────────────────────────────────────────────
    b_mask     = obs['Major_CT'] == 'CD4T'
    infer_mask = b_mask & (obs['day'] == '0') & (obs['year'] == '1')
    eval_mask  = b_mask & (obs['day'] == '0') & (obs['year'] == '2')

    print(f'  CD4T day=0 year=1: {infer_mask.sum()} cells', flush=True)
    print(f'  CD4T day=0 year=2: {eval_mask.sum()}  cells', flush=True)

    # ── select 10 shared donors ───────────────────────────────────────────────
    donors = select_donors(obs, infer_mask, eval_mask, n=N_DONORS)
    donor_map = {d: f'donor_{i}' for i, d in enumerate(sorted(donors))}
    print(f'  Selected donors: {donors}', flush=True)

    donor_mask = obs['donor_id'].isin(donors)
    infer_idx  = obs.index[infer_mask & donor_mask]
    eval_idx   = obs.index[eval_mask  & donor_mask]
    print(f'  Inference cells: {len(infer_idx)}', flush=True)
    print(f'  Evaluation cells: {len(eval_idx)}', flush=True)

    # ── load subsets to memory ────────────────────────────────────────────────
    print('Loading inference subset to memory...', flush=True)
    infer_raw = raw[infer_idx].to_memory()
    print('Loading evaluation subset to memory...', flush=True)
    eval_raw  = raw[eval_idx].to_memory()

    # ── downsample ────────────────────────────────────────────────────────────
    print(f'Downsampling to {CELLS_PER_DONOR} cells/donor...', flush=True)
    infer_raw = downsample_per_donor(infer_raw)
    eval_raw  = downsample_per_donor(eval_raw)
    print(f'  Inference after downsample: {infer_raw.n_obs}', flush=True)
    print(f'  Evaluation after downsample: {eval_raw.n_obs}', flush=True)

    # ── process ───────────────────────────────────────────────────────────────
    print('Processing inference...', flush=True)
    infer = process(infer_raw, donor_map)
    infer = add_metadata(infer, 'inference')
    print(infer, flush=True)

    print('Processing evaluation...', flush=True)
    evalu = process(eval_raw, donor_map)
    evalu = add_metadata(evalu, 'evaluation')
    print(evalu, flush=True)

    # ── write ─────────────────────────────────────────────────────────────────
    # ── rna_all (inference + evaluation concatenated) ─────────────────────────
    print('Building rna_all...', flush=True)
    import anndata as ad
    rna_all = ad.concat([infer, evalu], join='outer')
    rna_all.uns = infer.uns.copy()
    rna_all.uns['split'] = 'full'

    # ── write ─────────────────────────────────────────────────────────────────
    print(f'Writing {par["inference"]}...', flush=True)
    infer.write_h5ad(par['inference'])

    print(f'Writing {par["evaluation"]}...', flush=True)
    evalu.write_h5ad(par['evaluation'])

    print(f'Writing {par["rna_all"]}...', flush=True)
    rna_all.write_h5ad(par['rna_all'])

    print('\nDone.', flush=True)
    print(f'  Inference:  {infer.shape}')
    print(f'  Evaluation: {evalu.shape}')
    print(f'  RNA all:    {rna_all.shape}')


if __name__ == '__main__':
    main()
