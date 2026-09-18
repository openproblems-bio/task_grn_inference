#!/usr/bin/env python3
"""shoko — fully automatic coverage-layered GRN inference (geneRNIB component).

Thin viash wrapper around the bundled generator (resources/: gen_final_fast.py,
gen_cov_variant_v2.py, stream_core.py). No dataset-specific configuration:
  * expression layer: 'lognorm' when present, else 'X_norm' (auto-detected);
    an explicit --layer is honored only when that layer exists in the file
  * contrast routing: obs['is_control'] -> full method; derivable control
    labels (e.g. 'pre_*' perturbation labels) -> derived control -> full
    method; otherwise em-only reference arm (standardized cell-level
    |Pearson| network)
  * response granularity: per-donor (donor_env) when obs['donor_id'] exists

Full method: gene exclusions (near-constant + exact-duplicate columns),
donor-aware perturbation-response matrix, response-profile |Pearson| with
self-environment downdate, rank fusion with the cell-level |Pearson| network,
2M-edge candidate pool with pair-identity filter (|corr| < 0.99) and self-loop
exclusion, rule-R post-processing, TF-coverage layer (quota 10, layer cap
5,000, 50,000 edges total) and a retention (W) gate. Deterministic; the seed
argument is accepted for API compatibility but unused.
"""
import os
import sys

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/300BCG_rna.h5ad',
    'tf_all': 'resources_test/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/shoko_net.h5ad',
    'layer': 'lognorm',
    'max_n_links': 50000,
    'seed': 32,
    'dataset_id': 'op',
    'num_workers': 2,
    'temp_dir': 'output/temdir',
    'apply_tf_methods': True,
}
meta = {
    'name': 'shoko',
    'resources_dir': 'src/methods/shoko/resources',
}
## VIASH END

import h5py  # noqa: E402  (viash convention: imports may follow the VIASH par block)


def _resources_dir():
    """Locate the bundled generator modules (viash meta first, then script-relative)."""
    try:
        d = meta['resources_dir']
        if d and os.path.isfile(os.path.join(d, 'gen_final_fast.py')):
            return d
    except NameError:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (os.path.join(here, 'resources'), here):
        if os.path.isfile(os.path.join(d, 'gen_final_fast.py')):
            return d
    raise FileNotFoundError('shoko resources (gen_final_fast.py) not found')


def parse_args(par):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--rna', type=str)
    p.add_argument('--tf_all', type=str)
    p.add_argument('--prediction', type=str)
    p.add_argument('--layer', type=str)
    p.add_argument('--max_n_links', type=int)
    p.add_argument('--dataset_id', type=str)
    args, _ = p.parse_known_args()
    for k, v in vars(args).items():
        if v is not None:
            par[k] = v
    return par


def main(par):
    res = _resources_dir()
    sys.path.insert(0, res)
    # Bind the helper modules from the bundled resources BEFORE importing the
    # generator, so the generator's own fallback sys.path entries can never
    # shadow the bundled copies (sys.modules cache wins).
    import gen_cov_variant_v2  # noqa: F401
    import stream_core  # noqa: F401
    import gen_final_fast as gff

    argv = ['gen_final_fast.py',
            '--rna-override', par['rna'],
            '--tf-all', par['tf_all'],
            '--prediction', par['prediction']]

    with h5py.File(par['rna'], 'r') as f:
        avail = set(f['layers'].keys())
        try:
            raw = f['uns/dataset_id'][()]
            uns_did = raw.decode() if isinstance(raw, bytes) else str(raw)
        except Exception:
            uns_did = ''

    # layer: honor the requested layer when present, else let the generator
    # auto-detect ('lognorm' else 'X_norm')
    req_layer = par.get('layer') or ''
    if req_layer and req_layer in avail:
        argv += ['--layer', req_layer]
    elif req_layer:
        print(f'[shoko] requested layer {req_layer!r} absent from {sorted(avail)}; '
              f'falling back to auto-detection', flush=True)

    # dataset_id: an explicit non-default value wins; otherwise the h5ad's own
    # uns label is authoritative (same convention as the official methods)
    did = par.get('dataset_id') or ''
    if did and did != 'op':
        argv += ['--dataset-id', did]
    elif uns_did:
        argv += ['--dataset-id', uns_did]
    elif did:
        argv += ['--dataset-id', did]

    if int(par.get('max_n_links') or 50000) != 50000:
        print('[shoko] warning: max_n_links != 50000 is ignored; '
              'the shoko recipe emits a fixed 50,000-edge network', flush=True)

    sys.argv = argv
    gff.main()


if __name__ == '__main__':
    par = parse_args(par)
    main(par)
