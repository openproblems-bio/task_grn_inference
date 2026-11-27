import shutil
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

## VIASH START
par = {
    'dataset_id': 'op',
    'input_rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'input_atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
    'input_evaluation': 'resources_test/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'input_evaluation_de': None,
    'input_evaluation_sc': None,
    'input_ws_background': None,
    'metric_ids': None,
    'output_train': 'output/op/train.h5ad',
    'output_train_atac': 'output/op/train_atac.h5ad',
    'output_solution': 'output/op/solution.h5ad',
    'output_solution_de': None,
    'output_solution_sc': None,
    'output_ws_background': None,
}
## VIASH END

# Applicable metrics per dataset (copied from the legacy src/utils/config.py)
DATASETS_METRICS = {
    'replogle': ['regression', 'ws_distance', 'tf_recovery', 'tf_binding', 'sem', 'gs_recovery', 'vc'],
    'norman': ['regression', 'ws_distance', 'tf_binding', 'gs_recovery', 'vc', 'tf_recovery'],
    'nakatake': ['regression', 'gs_recovery', 'vc'],
    'op': ['regression', 'vc', 'replicate_consistency', 'tf_binding', 'sem', 'gs_recovery'],
    '300BCG': ['regression', 'vc', 'replicate_consistency', 'tf_binding', 'gs_recovery'],
    'parsebioscience': ['regression', 'vc', 'replicate_consistency', 'tf_binding', 'sem', 'gs_recovery'],
    'xaira_HEK293T': ['regression', 'ws_distance', 'tf_recovery', 'tf_binding', 'sem', 'gs_recovery', 'vc'],
    'xaira_HCT116': ['regression', 'ws_distance', 'tf_recovery', 'tf_binding', 'sem', 'gs_recovery', 'vc'],
    'MSCIC': ['regression', 'gs_recovery', 'replicate_consistency'],
    'soundlife': ['regression', 'replicate_consistency', 'tf_binding', 'gs_recovery'],
    'soundlife_vaccine': ['regression', 'replicate_consistency', 'tf_binding', 'gs_recovery'],
}
NORM_LAYERS = ('X_norm', 'lognorm')


def is_integer(x):
    d = x.data if sp.issparse(x) else np.asarray(x).ravel()
    return bool(np.allclose(d, np.round(d)))


def standardise(adata, dataset_id, with_counts=True, with_normalized=True):
    """Rename the dataset-specific normalized layer to `normalized`; counts are
    `layers['counts']` or X when X holds integers; otherwise counts are omitted."""
    layers = {}
    if with_normalized:
        src = next((l for l in NORM_LAYERS if l in adata.layers), None)
        if src is None:
            raise ValueError(f'None of {NORM_LAYERS} found in layers {list(adata.layers)}')
        norm = adata.layers[src]
        layers['normalized'] = norm.astype(np.float64)
    if with_counts:
        counts = adata.layers['counts'] if 'counts' in adata.layers else (adata.X if is_integer(adata.X) else None)
        if counts is not None:
            layers['counts'] = counts.astype(np.int32)
    uns = {k: v for k, v in adata.uns.items() if isinstance(v, str)}
    uns['dataset_id'] = dataset_id
    uns.setdefault('normalization_id', 'lognorm')
    if uns['normalization_id'] == 'logorm':  # typo in some source files
        uns['normalization_id'] = 'lognorm'
    return layers, uns


def build(adata, dataset_id, uns_extra=None, **kw):
    layers, uns = standardise(adata, dataset_id, **kw)
    uns.update(uns_extra or {})
    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    return ad.AnnData(obs=obs, var=adata.var.copy(), layers=layers, uns=uns)


dataset_id = par['dataset_id']
metric_ids = par['metric_ids'] or DATASETS_METRICS.get(dataset_id)
if not metric_ids:
    raise ValueError(f'No metric_ids for dataset {dataset_id}; pass --metric_ids')

print('>> train (RNA)', flush=True)
build(ad.read_h5ad(par['input_rna']), dataset_id).write_h5ad(par['output_train'], compression='gzip')

if par.get('input_atac') and par.get('output_train_atac'):
    print('>> train (ATAC)', flush=True)
    atac = ad.read_h5ad(par['input_atac'])
    build(atac, dataset_id, with_normalized=False).write_h5ad(par['output_train_atac'], compression='gzip')

print('>> solution (bulk)', flush=True)
build(ad.read_h5ad(par['input_evaluation']), dataset_id, with_counts=False,
      uns_extra={'metric_ids': list(metric_ids)}).write_h5ad(par['output_solution'], compression='gzip')

if par.get('input_evaluation_sc') and par.get('output_solution_sc'):
    print('>> solution (sc)', flush=True)
    build(ad.read_h5ad(par['input_evaluation_sc']), dataset_id, with_counts=False).write_h5ad(
        par['output_solution_sc'], compression='gzip')

if par.get('input_evaluation_de') and par.get('output_solution_de'):
    print('>> solution (DE)', flush=True)
    de = ad.read_h5ad(par['input_evaluation_de'])
    de.uns['dataset_id'] = dataset_id
    de.uns.setdefault('normalization_id', 'lognorm')
    de.write_h5ad(par['output_solution_de'], compression='gzip')

if par.get('input_ws_background') and par.get('output_ws_background'):
    print('>> ws background', flush=True)
    shutil.copy(par['input_ws_background'], par['output_ws_background'])
