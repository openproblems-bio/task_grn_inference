#!/usr/bin/env python3
"""eval23 gen_final_fast.py — simplified final-method generator (依頼23).

Final method (exact order, matching gen_warm.py + rule R pre-registration):
  1. W retention gate on the 50k ranking (build_variant output, pre-capfill):
     retention median > 0.4 -> W arm (per-source top-100 by s = weight *
     retention from the 50k net, e16 fill verbatim). W takes precedence.
  2. else rule R: candidate-pool em median >= 0.3 -> capfill(50k net) else
     the 50k net (BD) as-is.

Engineering (pre-registered, prereg_simplify.md):
  F2  streaming sufficient statistics when data is large (Mode B); exact
      dense path otherwise (Mode A). Numerically verified equivalent.
  F3  pool filter as em lookup (exact; removes the dominant compute stage).
Toggles (verification arms):
  --flat             S2: replace build_variant layers with flat top-50k
  --env-granularity  S3: force env-level response (pb donor_env bypass)

Usage: gen_final_fast.py <ds> [--flat] [--env-granularity] [--mode auto|A|B]
       gen_final_fast.py --rna-override X.h5ad --tf-all tf.csv --prediction out.h5ad
Output: /workspace/eval23/gen_fast/{ds}_{tag}.h5ad + _candidates.csv + _summary.json

eval25: outer gate (resolve_contrast) auto-routes em-only vs full method;
        ablation flags --no-exclusions / --no-pool-filter.
eval26: fully automatic — no per-dataset config. layer auto ('lognorm' else
        'X_norm'); granularity auto (donor_env when donor_id present);
        verified to reproduce all 12 benchmark networks bit-identically.
"""
import argparse
import hashlib
import json
import os
import sys
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, '/mnt/shared-workspace/xaira_checkpoint/code')
sys.path.insert(0, '/workspace/eval23')
import gen_cov_variant_v2 as gv  # noqa: E402
from stream_core import H5AdReader, StreamStats, streaming_randomized_pca  # noqa: E402

REPO = '/mnt/shared-workspace/task_grn_inference'
OUT = '/workspace/eval23/gen_fast'

# eval26: no per-dataset config — layer and granularity are auto-detected
# (layer: 'lognorm' else 'X_norm'; granularity: donor_env when donor_id is
# present). Verified to reproduce all 12 benchmark networks bit-identically.
MODE_A_MAX_BYTES = 8e9
N_PC = 5
W_GATE = 0.4
RULE_R = 0.3
CAP = 100
TOTAL = 50_000
TOP_N = 2_000_000
POOL_FACTOR = 3
PAIR_CORR_MAX = 0.99


def extract_topn_f3(fused, em_abs, tfs_present, tf_cols, genes, top_n,
                    pair_corr_max, pool_factor):
    """gv.extract_topn with the pair-correlation filter as an em lookup (F3).
    em_abs[i, j] = |Pearson(TF_i, gene_j)| on the same standardized cells."""
    n_tf, n_genes = fused.shape
    flat = fused.ravel()
    pool = min(int(top_n) * int(pool_factor), flat.size)
    idx = np.argpartition(-flat, pool - 1)[:pool]
    idx = idx[np.argsort(-flat[idx], kind='stable')]
    rows, cols = np.divmod(idx, n_genes)
    keep_mask = em_abs[rows, cols] < pair_corr_max
    keep_mask &= tf_cols[rows] != cols  # eval25: never emit self-loops (no-op when filter on)
    idx = idx[keep_mask][:top_n]
    rows, cols = np.divmod(idx, n_genes)
    tf_arr = np.asarray(tfs_present)
    return pd.DataFrame({'source': tf_arr[rows], 'target': genes[cols],
                         'weight': flat[idx]})


def capfill_e16(net, cand, cap=CAP):
    """Per-source top-cap by weight + fill m=1 per uncovered target (e16 rule,
    verbatim semantics as in gen_warm.py)."""
    capnet = (net.sort_values('weight', ascending=False, kind='stable')
                 .groupby('source', sort=True).head(cap).reset_index(drop=True))
    in_deg = capnet.groupby('target').size().to_dict()
    out_deg = capnet.groupby('source').size().to_dict()
    existing = set(zip(capnet['source'], capnet['target']))
    add = []
    for s_, t_, w_ in cand.sort_values('weight', ascending=False,
                                       kind='stable').itertuples(index=False):
        if in_deg.get(t_, 0) >= 1:
            continue
        if out_deg.get(s_, 0) >= cap or (s_, t_) in existing:
            continue
        add.append((s_, t_, w_))
        in_deg[t_] = 1
        out_deg[s_] = out_deg.get(s_, 0) + 1
    return pd.concat([capnet, pd.DataFrame(add, columns=['source', 'target', 'weight'])],
                     ignore_index=True)


def retention_stats(net, R, rC, rC2, tf_row, gcol, tf_cols):
    """gen_warm partial_stats: (r_ij, retention, ok) per network edge."""
    r = net['source'].map(tf_row).to_numpy()
    c = net['target'].map(gcol).to_numpy()
    ok = ~(pd.isna(r) | pd.isna(c))
    n = len(net)
    r_ij = np.full(n, np.nan, dtype=np.float32)
    ret = np.full(n, np.nan, dtype=np.float32)
    ri = r[ok].astype(int)
    ci = c[ok].astype(int)
    r_ij[ok] = R[ri, ci]
    gi = tf_cols[ri]
    num = r_ij[ok] - (rC[gi] * rC[ci]).sum(axis=1)
    den = np.sqrt((1 - rC2[gi]) * (1 - rC2[ci]))
    r_pc = num / den
    ret[ok] = np.clip(np.abs(r_pc) / np.maximum(np.abs(r_ij[ok]), 1e-9), 0, 1)
    return r_ij, ret, ok


CTRL_LABELS = {'control', 'ctrl', 'untreated', 'mock', 'vehicle', 'dmso',
               'baseline', 'healthy', 'normal'}


def resolve_contrast(obs, gran):
    """eval25 outer gate (PLAN.md 依頼25). obs: pandas DataFrame.
    Returns (labels, ctrl, donors, gate). gate['route']:
      full    — is_control present and non-degenerate (unchanged behavior)
      derived — is_control missing; control-like labels found in perturbation
                (CTRL_LABELS case-insensitive, or 'pre_'/'pre-'/'pre'); requires
                both groups non-empty and (donor_env) per-donor coverage
      em_only — no usable intervention contrast
    """
    gate = {'route': None, 'reason': '', 'control_labels': [],
            'n_ctrl': 0, 'n_treated': 0, 'donor_coverage_ok': None}
    donors = obs['donor_id'].astype(str).to_numpy() if 'donor_id' in obs.columns else None
    if 'perturbation' not in obs.columns:
        gate.update(route='em_only', reason='no perturbation column')
        return None, None, donors, gate
    labels = obs['perturbation'].astype(str).to_numpy()
    if 'is_control' in obs.columns:
        ctrl = obs['is_control'].astype(bool).to_numpy()
        if ctrl.sum() > 0 and (~ctrl).sum() > 0 and len(set(labels[~ctrl])) >= 1:
            gate.update(route='full', n_ctrl=int(ctrl.sum()),
                        n_treated=int((~ctrl).sum()))
            return labels, ctrl, donors, gate
        gate.update(route='em_only', reason='degenerate is_control')
        return None, None, donors, gate
    # aggressive derivation (user-approved): control-like label patterns
    lab_lower = [lab.lower() for lab in labels]
    is_ctrl_like = np.array([lab in CTRL_LABELS or lab == 'pre'
                             or lab.startswith('pre_') or lab.startswith('pre-')
                             for lab in lab_lower])
    if is_ctrl_like.sum() == 0 or (~is_ctrl_like).sum() == 0:
        gate.update(route='em_only',
                    reason='no control-like or no treated-like labels')
        return None, None, donors, gate
    if gran == 'donor_env' and donors is not None:
        cov = [bool(((donors == d) & is_ctrl_like).sum() > 0
                    and ((donors == d) & ~is_ctrl_like).sum() > 0)
               for d in set(donors)]
        gate['donor_coverage_ok'] = all(cov)
        if not gate['donor_coverage_ok']:
            gate.update(route='em_only',
                        reason='donor coverage incomplete for derived control')
            return None, None, donors, gate
    gate.update(route='derived',
                control_labels=sorted(set(labels[is_ctrl_like])),
                n_ctrl=int(is_ctrl_like.sum()),
                n_treated=int((~is_ctrl_like).sum()))
    return labels, is_ctrl_like, donors, gate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ds', nargs='?', default=None,
                    help='dataset id for repo-layout I/O; optional with --rna-override/--prediction')
    ap.add_argument('--flat', action='store_true')
    ap.add_argument('--env-granularity', action='store_true')
    ap.add_argument('--em-only', action='store_true',
                    help='eval24 Phase B reference: expression-correlation-only '
                         'degenerate variant (no response channel, no gates)')
    ap.add_argument('--no-exclusions', action='store_true',
                    help='eval25 ablation: disable std-floor + duplicate-gene exclusions')
    ap.add_argument('--no-pool-filter', action='store_true',
                    help='eval25 ablation: disable the pair-corr candidate-pool filter')
    ap.add_argument('--mode', default='auto', choices=['auto', 'A', 'B'])
    ap.add_argument('--rna-override', default=None,
                    help='optional alternate rna h5ad path (e.g. synthetic stress-test data)')
    ap.add_argument('--layer', default=None,
                    help="override expression layer (default: auto — 'lognorm' else 'X_norm')")
    ap.add_argument('--granularity', default=None, choices=['env', 'donor_env'],
                    help='override response granularity (default: auto — donor_env when donor_id present)')
    ap.add_argument('--tf-all', default=None, help='override tf_all.csv path')
    ap.add_argument('--prediction', default=None,
                    help='override output prediction h5ad path (component mode)')
    ap.add_argument('--dataset-id', default=None, help='override uns dataset_id')
    args = ap.parse_args()
    ds = args.ds or args.dataset_id or 'custom'
    t_all = time.time()
    stages = {}

    def mark(name, t0):
        stages[name] = round(time.time() - t0, 2)
        print(f'[{ds}/{tag}] {name}: {stages[name]}s', flush=True)

    # ---------------------------------------------------------- load -------
    t0 = time.time()
    rna_path = args.rna_override or f'{REPO}/resources/grn_benchmark/inference_data/{ds}_rna.h5ad'
    # shape-only probe via h5py (anndata backed='r' still loads LAYERS eagerly)
    with h5py.File(rna_path, 'r') as f:
        avail_layers = list(f['layers'].keys())
        layer = args.layer or ('lognorm' if 'lognorm' in avail_layers
                               else 'X_norm' if 'X_norm' in avail_layers else None)
        if layer is None:
            raise ValueError(f'no usable layer in {rna_path} (have {avail_layers})')
        try:
            raw = f['uns/dataset_id'][()]
            uns_dataset_id = raw.decode() if isinstance(raw, bytes) else str(raw)
        except Exception:
            uns_dataset_id = ''
        node = f[f'layers/{layer}']
        if isinstance(node, h5py.Group):
            n, G_all = (int(v) for v in node.attrs['shape'])
        else:
            n, G_all = int(node.shape[0]), int(node.shape[1])
    est_bytes = n * G_all * 4
    mode = args.mode if args.mode != 'auto' else ('A' if est_bytes < MODE_A_MAX_BYTES else 'B')
    print(f'[{ds}] mode={mode} cells={n} genes={G_all} layer={layer}', flush=True)

    # ---- load obs + eval25 outer gate (auto em-only/full routing) ----
    if mode == 'A':
        adata = ad.read_h5ad(rna_path)
        genes = adata.var_names.astype(str).to_numpy()
        obs_df = adata.obs
    else:
        reader = H5AdReader(rna_path, layer)
        genes = reader.genes
        obs_df = pd.DataFrame({c: reader.obs_col(c)
                               for c in ('perturbation', 'is_control', 'donor_id', 'cell_count')
                               if reader.has_obs(c)})
    # eval26: granularity auto — donor_env when donor_id exists (contrast route)
    gran = (args.granularity
            or ('env' if args.env_granularity
                else ('donor_env' if 'donor_id' in obs_df.columns else 'env')))
    if args.em_only:
        assert mode == 'A', 'em-only reference arm supports Mode A only'
        labels = ctrl = donors = weights = None
        gate = {'route': 'forced_em_only'}
    else:
        labels, ctrl, donors, gate = resolve_contrast(obs_df, gran)
        weights = None
        if gate['route'] != 'em_only':
            weights = (obs_df['cell_count'].to_numpy(dtype=np.float64)
                       if 'cell_count' in obs_df.columns
                       else np.ones(len(obs_df), dtype=np.float64))
    em_only = bool(args.em_only or gate['route'] == 'em_only')
    assert not (em_only and mode == 'B'), 'em-only (auto or forced) supports Mode A only'
    tag = ('emonly' if em_only else
           'flat' if args.flat else ('envcap' if args.env_granularity else 'fast'))
    if args.flat and em_only:
        tag += '_flat'
    if args.no_exclusions:
        tag += '_noex'
    if args.no_pool_filter:
        tag += '_nofilt'
    if args.mode == 'B':
        tag += '_modeB'
    print(f'[{ds}/{tag}] outer gate: {gate}', flush=True)
    mark('load_meta', t0)

    tf_all = pd.read_csv(args.tf_all or f'{REPO}/resources/grn_benchmark/prior/tf_all.csv',
                         header=None)[0].astype(str)
    tf_set = set(tf_all)
    tf_cols_full = np.array([i for i, g in enumerate(genes) if g in tf_set])

    if mode == 'A':
        t0 = time.time()
        mat = adata.layers[layer]
        if hasattr(mat, 'toarray'):
            mat = mat.toarray()
        x = np.asarray(mat, dtype=np.float32)
        del mat
        del adata
        mark('densify', t0)

        t0 = time.time()
        if args.no_exclusions:
            floor, dup = set(), set()
        else:
            floor, dup = gv.find_exclusions(x, genes, 0.0032)
        mark('exclusions', t0)

        t0 = time.time()
        if em_only:
            r = None
        elif gran == 'donor_env':
            r, envs, n_ctrl, row_donors = gv.build_donor_env_response_fast(x, labels, donors, ctrl, weights)
        else:
            r, envs, n_ctrl = gv.build_env_response(x, labels, ctrl)
        mark('response', t0)

        keep = np.array([g not in (floor | dup) for g in genes])
        genes_k = genes[keep]
        t0 = time.time()
        z = x[:, keep].astype(np.float32)
        z -= z.mean(axis=0)
        sd = z.std(axis=0)
        sd[sd == 0] = 1.0
        z /= sd
        mark('zscore', t0)
        del x
    else:
        # -------------------------------------------------- Mode B: streaming
        t0 = time.time()
        if gran == 'donor_env':
            group_keys = np.array([f'{d}|{e}' for d, e in zip(donors, labels)])
            ss = StreamStats(G_all, tf_cols_full, group_keys, ctrl,
                             donors=donors, weights=weights)
        else:
            ss = StreamStats(G_all, tf_cols_full, labels, ctrl)
        hashes = [hashlib.md5() for _ in range(G_all)]
        CH = 10_000
        for s in range(0, n, CH):
            e = min(s + CH, n)
            chunk = reader.layer_rows(s, e)  # dense float32 straight from disk
            ss.update(chunk, np.arange(s, e))
            for gi in range(G_all):
                hashes[gi].update(np.ascontiguousarray(chunk[:, gi]).tobytes())
        mark('stream_pass1', t0)

        t0 = time.time()
        mean, std = ss.mean_std()
        floor = set(genes[std <= 0.0032])
        seen, dup = {}, set()
        for gi, g in enumerate(genes):
            h = hashes[gi].hexdigest()
            if h in seen:
                dup.add(g)
            else:
                seen[h] = g
        dup -= floor
        mark('exclusions', t0)

        keep = np.array([g not in (floor | dup) for g in genes])
        genes_k = genes[keep]
        r = ss.response_matrix()
        envs = [k.split('|', 1)[1] if gran == 'donor_env' else k for k in ss.g2i.keys()]
        if gran == 'donor_env':
            # envs must match build_donor_env_response_fast ordering: sorted (donor, env)
            order = np.argsort(list(ss.g2i.keys()))
            r = r[order]
            envs = [list(ss.g2i.keys())[i].split('|', 1)[1] for i in order]
        mark('response', t0)

    # --------------------------------------------------- shared downstream -
    t0 = time.time()
    keep_idx = np.flatnonzero(keep)
    tfs_present = [g for g in genes_k if g in tf_set]
    tf_cols = np.array([i for i, g in enumerate(genes_k) if g in tf_set])
    if em_only:
        sm = None
        mark('score_matrix_skipped_em_only', t0)
    else:
        assert r is not None  # response matrix is always built on the full-method path
        y = r[:, keep]
        y = y - y.mean(axis=0)
        env_row = {}
        for i, e2 in enumerate(envs):
            env_row.setdefault(e2, []).append(i)
        env_set = set(envs)
        self_env_rows = {i: env_row[t] for i, t in enumerate(tfs_present) if t in env_set}
        sm = gv.score_matrix(y, tf_cols, self_env_rows)
        mark('score_matrix', t0)
        del y, r

    t0 = time.time()
    if mode == 'A':
        em_signed = None
        em = gv.expr_corr_matrix(z, tf_cols)  # |r|
    else:
        em_full = ss.em_signed()              # TF x G_all signed
        tf_row_full = {g: i for i, g in enumerate(genes[tf_cols_full])}
        rows_keep = np.array([tf_row_full[g] for g in tfs_present])
        em_signed = em_full[np.ix_(rows_keep, keep_idx)]
        em = np.abs(em_signed)
    mark('em', t0)

    t0 = time.time()
    if em_only:
        fused = em.copy()  # degenerate reference: em ranking is the only channel
    else:
        fused = gv.rank_fuse(sm, em)
    mark('rank_fuse', t0)
    del sm

    t0 = time.time()
    pc_max = 1.0 if args.no_pool_filter else PAIR_CORR_MAX
    cand = extract_topn_f3(fused, em, tfs_present, tf_cols, genes_k,
                           TOP_N, pc_max, POOL_FACTOR)
    mark('extract_topn_F3', t0)

    # rule R: candidate-pool em median
    tf_row = {g: i for i, g in enumerate(tfs_present)}
    gcol = {g: i for i, g in enumerate(genes_k)}
    ri = cand['source'].map(tf_row).to_numpy()
    ci = cand['target'].map(gcol).to_numpy()
    cand_em_p50 = float(np.median(em[ri, ci]))
    if em_only:
        post = 'em_only'  # gates are not applicable to the degenerate reference
    else:
        post = 'capfill' if cand_em_p50 >= RULE_R else 'BD'
    print(f'[{ds}/{tag}] rule R: cand em p50={cand_em_p50:.4f} -> {post}', flush=True)

    t0 = time.time()
    if args.flat:
        net = cand.head(TOTAL).reset_index(drop=True)
        info = {'L_sat': None, 'L': None, 'core_edges': TOTAL}
    else:
        base = cand.head(TOTAL)
        net, info = gv.build_variant(base, cand, 10, 5000, TOTAL, layer_weight='real')
    mark('build_variant', t0)

    # ----------------------- W retention gate (on the 50k net, pre-capfill) -
    t0 = time.time()
    if mode == 'A':
        from sklearn.utils.extmath import randomized_svd
        U, S, Vt = randomized_svd(z, n_components=N_PC, random_state=0)
    else:
        mean, std = ss.mean_std()
        std_k = std[keep]
        mean_k = mean[keep]

        def load_chunk(a, b):
            return reader.layer_rows(a, b)[:, keep_idx]
        U, S, Vt = streaming_randomized_pca(load_chunk, n, len(genes_k),
                                            mean_k, std_k, n_components=N_PC)
    P = U * S
    pstd = P.std(axis=0)
    if mode == 'A':
        rC = ((z.T @ P) / n / pstd).astype(np.float32)
        zt = np.ascontiguousarray(z[:, tf_cols].T)
        R = ((zt @ z) / n).astype(np.float32)
        del zt
    else:
        # rC = corr(gene, PC): (Z^T P)/n/pstd via chunked accumulation
        rC = np.zeros((len(genes_k), N_PC), dtype=np.float64)
        std_safe_k = np.where(std_k == 0, 1.0, std_k)
        for s in range(0, n, 10_000):
            e = min(s + 10_000, n)
            zc = load_chunk(s, e).astype(np.float64)
            zc -= mean_k
            zc /= std_safe_k
            rC += zc.T @ P[s:e]
        rC = (rC / n / pstd).astype(np.float32)
        R = em_signed
    rC2 = np.clip((rC ** 2).sum(axis=1), 0, 1 - 1e-6)
    r_ij, ret, ok = retention_stats(net, R, rC, rC2, tf_row, gcol, tf_cols)
    ret_med = float(np.nanmedian(ret))
    mark('w_gate', t0)

    w_applied = False
    if em_only:
        pass  # reference arm: W gate computed as diagnostic but never applied
    elif ret_med > W_GATE:
        t0 = time.time()
        ret_use = np.where(ok, ret, 1.0)
        net2 = net.copy()
        net2['s'] = net2['weight'].to_numpy() * ret_use
        capnet = (net2.sort_values('s', ascending=False, kind='stable')
                      .groupby('source', sort=True).head(CAP)
                      .reset_index(drop=True)[['source', 'target', 'weight']])
        net = capfill_e16(capnet, cand)
        w_applied = True
        mark('w_reselect', t0)
    elif post == 'capfill':
        t0 = time.time()
        net = capfill_e16(net, cand)
        mark('post_capfill', t0)
    print(f'[{ds}/{tag}] W gate: retention median={ret_med:.4f} -> W={w_applied}', flush=True)

    # ------------------------------------------------------------- write ---
    out = ad.AnnData(np.zeros((0, 0), dtype=np.float32))
    out.uns['dataset_id'] = (args.dataset_id or (ds if ds != 'custom' else '')
                             or uns_dataset_id or 'unknown')
    out.uns['method_id'] = 'shoko'
    out.uns['prediction'] = net.assign(weight=net['weight'].astype(str))
    path = args.prediction or f'{OUT}/{ds}_{tag}.h5ad'
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    out.write_h5ad(path)
    base = path[:-len('.h5ad')] if path.endswith('.h5ad') else path
    cand.to_csv(f'{base}_candidates.csv', index=False)
    summary = {'dataset': ds, 'tag': tag, 'mode': mode, 'layer': layer,
               'granularity': gran,
               'flat': args.flat, 'em_only': em_only, 'outer_gate': gate,
               'no_exclusions': args.no_exclusions,
               'no_pool_filter': args.no_pool_filter,
               'cand_em_p50': cand_em_p50, 'post': post,
               'retention_median': ret_med, 'w_applied': w_applied,
               'final_edges': len(net), 'stages': stages,
               'total_sec': round(time.time() - t_all, 2)}
    with open(f'{base}_summary.json', 'w') as f:
        json.dump(summary, f, indent=1)
    print(f'[{ds}/{tag}] DONE edges={len(net)} total={summary["total_sec"]}s', flush=True)


if __name__ == '__main__':
    main()
