"""
SoundLife sc B-cell processing — Scenario 2.

Temporal split: GRN inferred from pre- and early post-vaccination states
(year 1 only), evaluated on the 90-day memory consolidation timepoint
across year 1 and year 2.

  Inference:  B cells, day=0 + day=7, year=1          (10 donors)
  Evaluation: B cells, day=90,        year=1 + year=2 (same 10 donors)

Goal: evaluate whether a GRN inferred from baseline and acute vaccination
response (day 7, year 1) correctly predicts transcriptional state at day 90
(memory B cell formation / immune consolidation) in year 1 and year 2.

Source: /vol/projects/jnourisa/hiara/datasets_old/sc/soundlife.h5ad
Outputs:
  resources/grn_benchmark/inference_data/soundlife_B_rna_s2.h5ad
  resources/grn_benchmark/evaluation_data/soundlife_B_sc_s2.h5ad
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

RAW      = '/vol/projects/jnourisa/hiara/datasets_old/sc/soundlife.h5ad'
N_DONORS = 10

par = {
    'inference':  os.path.join(ROOT, 'resources/grn_benchmark/inference_data/soundlife_vaccine_rna.h5ad'),
    'evaluation': os.path.join(ROOT, 'resources/grn_benchmark/evaluation_data/soundlife_vaccine_sc.h5ad'),
    'rna_all':    os.path.join(ROOT, 'resources/extended_data/soundlife_vaccine_rna_all.h5ad'),
}

PERTURBATION_MAP = {
    '0':  'pre_vaccination',
    '7':  'post_vaccination_d7',
    '90': 'post_vaccination_d90',
}


def select_donors(obs, infer_mask, eval_mask, n=N_DONORS, seed=0):
    """Return n donor_ids present in both inference and evaluation."""
    infer_donors = set(obs.loc[infer_mask, 'donor_id'].unique())
    eval_donors  = set(obs.loc[eval_mask,  'donor_id'].unique())
    shared = sorted(infer_donors & eval_donors)
    rng = np.random.default_rng(seed)
    return list(rng.choice(shared, size=min(n, len(shared)), replace=False))


def harmonise_obs(obs, donor_map):
    out = pd.DataFrame(index=obs.index)
    out['donor_id']       = obs['donor_id'].map(donor_map)
    out['cell_type']      = 'B'
    out['cell_type_orig'] = obs['Sub_CT'].astype(str)
    out['perturbation']   = obs['day'].astype(str).map(PERTURBATION_MAP).fillna('unknown')
    out['day']            = obs['day'].astype(str)
    out['year']           = obs['year'].astype(str)
    out['visitName']      = obs['visitName'].astype(str)
    out['donor_age']      = obs['sample.subjectAgeAtDraw'].astype(str)
    out['donor_sex']      = obs['subject.biologicalSex'].astype(str)
    out.index.name = 'obs_id'
    return out


def add_metadata(adata, split_label):
    adata.uns['dataset_id']       = 'soundlife'
    adata.uns['dataset_name']     = 'SoundLife'
    adata.uns['dataset_summary']  = (
        'SoundLife longitudinal flu-vaccination cohort. '
        'PBMCs from healthy adults sampled at day 0 (pre-vax), day 7 (acute response), '
        'and day 90 (memory consolidation) across two flu seasons.'
    )
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    adata.uns['split']            = split_label
    adata.uns['scenario']         = '2'
    adata.uns['scenario_description'] = (
        'Inference: B cells, day=0+7, year=1. '
        'Evaluation: B cells, day=90, year=1+2 (same donors).'
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
    b_mask     = obs['Major_CT'] == 'B'
    year_mask  = obs['year'].isin(['1', '2'])
    infer_mask = b_mask & (obs['year'] == '1') & obs['day'].isin(['0', '7'])
    eval_mask  = b_mask & year_mask & (obs['day'] == '90')

    print(f'  B day=0+7 year=1:   {infer_mask.sum()} cells', flush=True)
    print(f'  B day=90  year=1+2: {eval_mask.sum()}  cells', flush=True)

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

    # ── process ───────────────────────────────────────────────────────────────
    print('Processing inference...', flush=True)
    infer = process(infer_raw, donor_map)
    infer = add_metadata(infer, 'inference')
    print(infer, flush=True)

    print('Processing evaluation...', flush=True)
    evalu = process(eval_raw, donor_map)
    evalu = add_metadata(evalu, 'evaluation')
    print(evalu, flush=True)

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
