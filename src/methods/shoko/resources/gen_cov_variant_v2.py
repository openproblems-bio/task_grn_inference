#!/usr/bin/env python3
"""gen_cov_variant_v2.py — covlayer improvement-round generator (Request 14).

Extends the bit-validated covlayer recipe (component_covlayer_script.py) with
four orthogonal axes, each motivated by a Request-13 diagnostic:
  --granularity {env,donor_env}   A: donor x env pseudobulk response points
                                  (per-donor control subtraction, cell_count-
                                  weighted well means). Falls back to env when
                                  obs has no donor_id.
  --signal {pairwise,ridge}       C: per-target ridge of target response on all
                                  TF responses (standardized coefficients,
                                  alpha by GCV over a fixed grid, exact reduced
                                  solve for TF targets excluding own column).
  --layer-weight {artificial,real} B: layer edges keep their real candidate
                                  scores instead of the artificial decay ramp.
  --fuse-expr {none,rank50}       D: 50/50 rank fusion of the signal matrix
                                  with |cell-level expression Pearson corr|.

Default mode (env, pairwise, artificial, none) is byte-identical in logic to
component_covlayer_script.py (Gate 1 verifies against C0_covLstar_ranking.csv).

Outputs: <out>_ranking.csv (50k variant), <out>_candidates.csv (top-n candidate
ranking), <out>_summary.json, and a prediction h5ad (method_id='shoko').

Usage:
  python gen_cov_variant_v2.py --dataset DS --rna H5AD --tf_all TF_CSV \
      --granularity env --signal pairwise --layer-weight artificial \
      --fuse-expr none --out-prefix OUT --prediction PRED.h5ad
"""
import argparse
import json
import time

import anndata as ad
import numpy as np
import pandas as pd

RIDGE_ALPHAS = np.logspace(-3, 3, 13)  # fixed GCV grid (not tuned)


# ---------------------------------------------------------------- response --
def build_env_response(x, labels, ctrl_mask):
    """Original env-level response (unweighted well means, global control)."""
    if ctrl_mask.sum() == 0:
        raise ValueError('No control cells (obs["is_control"] all False)')
    ctrl_mean = x[ctrl_mask].mean(axis=0)
    envs, rows = [], []
    for env in pd.unique(labels[~ctrl_mask]):
        idx = np.flatnonzero(labels == env)
        envs.append(str(env))
        rows.append(x[idx].mean(axis=0) - ctrl_mean)
    return np.vstack(rows).astype(np.float64), envs, int(ctrl_mask.sum())


def build_donor_env_response_fast(x, labels, donors, ctrl_mask, weights):
    """A: donor x env pseudobulk response, per-donor control subtraction.

    Point(d, e) = cell_count-weighted mean of wells in (d, e) minus the
    cell_count-weighted mean of donor d's control wells. Rows ordered by
    sorted (donor, env); env labels kept per row for self-env exclusion.
    """
    donors = np.asarray(donors)
    labels = np.asarray(labels)
    w = weights.astype(np.float64)
    xw = x.astype(np.float64) * w[:, None]
    # per-donor control mean
    ctrl_mean_d = {}
    for d in np.unique(donors):
        m = (donors == d) & ctrl_mask
        ctrl_mean_d[d] = xw[m].sum(axis=0) / w[m].sum()
    combos = sorted(set(zip(donors[~ctrl_mask], labels[~ctrl_mask])))
    envs, rows, row_donors = [], [], []
    for d, e in combos:
        m = (donors == d) & (labels == e)
        rows.append(xw[m].sum(axis=0) / w[m].sum() - ctrl_mean_d[d])
        envs.append(str(e))
        row_donors.append(d)
    return (np.vstack(rows).astype(np.float64), envs,
            int(ctrl_mask.sum()), row_donors)


# --------------------------------------------------------------- exclusions --
def find_exclusions(x, genes, floor_std):
    import hashlib
    std = x.std(axis=0)
    floor = set(genes[std <= floor_std])
    seen, dup_excess = {}, set()
    for i, g in enumerate(genes):
        h = hashlib.md5(np.ascontiguousarray(x[:, i]).tobytes()).hexdigest()
        if h in seen:
            dup_excess.add(g)
        else:
            seen[h] = g
    dup_excess -= floor
    return floor, dup_excess


# ----------------------------------------------------------------- signals --
def score_matrix(y, tf_cols, self_env_rows):
    """Original pairwise |Pearson| of response profiles (exact copy)."""
    y_tf = y[:, tf_cols]
    cov = y_tf.T @ y
    var_t = (y_tf ** 2).sum(axis=0)
    var_g = (y ** 2).sum(axis=0)
    scores = np.abs(cov) / np.sqrt(var_t[:, None] * var_g[None, :])
    for row, env_rows in self_env_rows.items():
        r = y[env_rows, :]
        cov_r = cov[row] - (y_tf[env_rows, row][:, None] * r).sum(axis=0)
        var_t_t = var_t[row] - float((y_tf[env_rows, row] ** 2).sum())
        var_g_g = var_g - (r ** 2).sum(axis=0)
        scores[row] = np.abs(cov_r) / np.sqrt(var_t_t * var_g_g)
    return scores


def ridge_scores(y, tf_cols, alphas=RIDGE_ALPHAS):
    """C: per-target ridge of standardized target response on TF responses.

    Returns (score matrix |beta|, per-target chosen alpha). y must be centered
    per gene. TF targets are fit without their own predictor column (exact
    reduced Cholesky solve at the GCV-selected alpha).
    """
    n, p_all = y.shape
    sd = y.std(axis=0)
    sd[sd == 0] = 1.0
    ys = y / sd
    X = ys[:, tf_cols]
    p = X.shape[1]
    G = X.T @ X
    P = X.T @ ys                       # p x n_genes
    ynorm = (ys ** 2).sum(axis=0)      # per target
    lam, Q = np.linalg.eigh(G)
    QP = Q.T @ P                       # p x n_genes
    tr_H = np.array([np.sum(lam / (lam + a)) for a in alphas])
    # GCV per alpha per target via sufficient statistics
    gcv = np.empty((len(alphas), p_all))
    for ai, a in enumerate(alphas):
        d = QP / (lam + a)[:, None]                    # p x n_genes
        b_p = (QP * d).sum(axis=0)                     # beta^T p
        b_G_b = (lam[:, None] * (QP ** 2) / (lam + a)[:, None] ** 2).sum(axis=0)
        rss = ynorm - 2 * b_p + b_G_b
        rss = np.maximum(rss, 0.0)
        gcv[ai] = n * rss / (n - tr_H[ai]) ** 2
    best = np.argmin(gcv, axis=0)
    alpha_sel = alphas[best]
    # beta at selected alpha (full predictor set)
    d_sel = QP / (lam[:, None] + alpha_sel[None, :])   # (p, n_genes)
    B = Q @ d_sel                                      # p x n_genes
    # exact reduced solve for TF targets (drop own column)
    tf_set = {int(c): i for i, c in enumerate(tf_cols)}
    for j, i in tf_set.items():
        a = alpha_sel[j]
        keep = np.arange(p) != i
        Gk = G[np.ix_(keep, keep)]
        pk = P[keep, j]
        B[keep, j] = np.linalg.solve(Gk + a * np.eye(p - 1), pk)
        B[i, j] = 0.0  # self coefficient removed
    scores = np.abs(B)
    return scores, alpha_sel


def expr_corr_matrix(z, tf_cols):
    """|Pearson| of cell-level standardized expression, TF x genes."""
    n = z.shape[0]
    return np.abs((z[:, tf_cols].T @ z) / n)


def rank_fuse(sm, em):
    """D: 50/50 average of normalized ranks over the full TF x gene matrix."""
    flat_s, flat_e = sm.ravel(), em.ravel()
    N = flat_s.size
    r_s = np.empty(N, dtype=np.float64)
    r_s[np.argsort(flat_s, kind='stable')] = np.arange(N)
    r_e = np.empty(N, dtype=np.float64)
    r_e[np.argsort(flat_e, kind='stable')] = np.arange(N)
    fused = 0.5 * (r_s + r_e) / N
    return fused.reshape(sm.shape)


# ------------------------------------------------------------- extraction --
def extract_topn(sm, tfs_present, tf_cols, genes, top_n, z_cells, pair_corr_max, pool_factor):
    n_tf, n_genes = sm.shape
    flat = sm.ravel()
    pool = min(int(top_n) * int(pool_factor), flat.size)
    idx = np.argpartition(-flat, pool - 1)[:pool]
    idx = idx[np.argsort(-flat[idx], kind='stable')]
    rows, cols = np.divmod(idx, n_genes)
    src_gene_cols = tf_cols[rows]
    n_cells = z_cells.shape[0]
    keep_mask = np.ones(len(idx), dtype=bool)
    chunk = 4096
    for s in range(0, len(idx), chunk):
        a = z_cells[:, src_gene_cols[s:s + chunk]]
        b = z_cells[:, cols[s:s + chunk]]
        c = np.einsum('ij,ij->j', a, b) / n_cells
        keep_mask[s:s + chunk] = np.abs(c) < pair_corr_max
    idx = idx[keep_mask][:top_n]
    rows, cols = np.divmod(idx, n_genes)
    tf_arr = np.asarray(tfs_present)
    return pd.DataFrame({'source': tf_arr[rows], 'target': genes[cols], 'weight': flat[idx]})


# ------------------------------------------------------------------ layer --
def l_sat(ranking_sources, base, core_n, quota):
    avail = ranking_sources.value_counts()
    cc = base.head(core_n).groupby('source').size().reindex(avail.index).fillna(0)
    need = (quota - cc).clip(lower=0)
    need = need.where(avail >= quota, (avail - cc).clip(lower=0))
    return int(need.sum())


def l_sat_dynamic(ranking_sources, base, total, quota):
    L = l_sat(ranking_sources, base, min(45000, total), quota)
    for _ in range(50):
        L_new = l_sat(ranking_sources, base, total - L, quota)
        if L_new == L:
            break
        L = L_new
    return L


def build_variant(base, cand, quota, layer_cap, total, layer_weight='artificial'):
    L_sat = l_sat_dynamic(cand['source'], base, total, quota)
    L = min(int(layer_cap), L_sat)
    core_n = total - L
    core = base.head(core_n)
    core_pairs = set(zip(core['source'], core['target']))

    in_core = pd.MultiIndex.from_arrays([cand['source'], cand['target']]).isin(core_pairs)
    pool = cand[~in_core]
    pool_rank = np.flatnonzero(~in_core)

    cc = core.groupby('source').size()
    prio = pd.Series(pool_rank, index=pool['source'].values).groupby(level=0).min()
    need = (quota - cc).clip(lower=0)
    eligible = prio[need.reindex(prio.index).fillna(quota) > 0].sort_values()

    pool_by_src = dict(tuple(pool.groupby('source', sort=False)))
    rows, budget = [], L
    for s in eligible.index:
        if budget <= 0:
            break
        n = int(min(need.get(s, quota), budget))
        take = pool_by_src[s].head(n)
        rows.append(take[['source', 'target', 'weight']])
        budget -= len(take)
    lay = (pd.concat(rows).reset_index(drop=True) if rows
           else pd.DataFrame(columns=['source', 'target', 'weight']))

    core_min = float(core['weight'].astype(float).abs().min())
    if layer_weight == 'artificial':
        lay = lay.assign(weight=core_min * (0.5 - 0.25 * np.arange(len(lay)) / max(len(lay), 1)))
    else:  # real: keep candidate scores (all <= core_min by construction)
        lay = lay.assign(weight=lay['weight'].astype(float))

    variant = (pd.concat([core, lay], ignore_index=True)
                 .sort_values('weight', ascending=False)
                 .reset_index(drop=True))
    info = {'L_sat': int(L_sat), 'L': int(L), 'core_edges': int(len(core)),
            'core_min_abs_weight': core_min, 'sources': int(variant['source'].nunique()),
            'layer_weight_mode': layer_weight,
            'layer_weight_min': float(lay['weight'].min()) if len(lay) else None,
            'layer_weight_max': float(lay['weight'].max()) if len(lay) else None}
    return variant, info


# ------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', required=True)
    ap.add_argument('--rna', required=True)
    ap.add_argument('--tf_all', required=True)
    ap.add_argument('--prediction', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--layer', default='lognorm')
    ap.add_argument('--max-n-links', type=int, default=50000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--quota', type=int, default=10)
    ap.add_argument('--layer-cap', type=int, default=5000)
    ap.add_argument('--floor-std', type=float, default=0.0032)
    ap.add_argument('--pair-corr-max', type=float, default=0.99)
    ap.add_argument('--pool-factor', type=int, default=3)
    ap.add_argument('--top-n-candidates', type=int, default=2000000)
    ap.add_argument('--granularity', choices=['env', 'donor_env'], default='env')
    ap.add_argument('--signal', choices=['pairwise', 'ridge'], default='pairwise')
    ap.add_argument('--layer-weight', choices=['artificial', 'real'], default='artificial')
    ap.add_argument('--fuse-expr', choices=['none', 'rank50'], default='none')
    ap.add_argument('--tf-out-cap', type=int, default=0,
                    help='if >0, keep only top-K targets per TF (by weight) in the final ranking')
    args = ap.parse_args()
    t0 = time.time()

    adata = ad.read_h5ad(args.rna)
    if args.layer not in adata.layers:
        raise ValueError(f'Layer {args.layer!r} not found')
    if 'is_control' not in adata.obs or 'perturbation' not in adata.obs:
        raise ValueError('requires obs["is_control"] and obs["perturbation"]')

    mat = adata.layers[args.layer]
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    x = np.asarray(mat, dtype=np.float32)
    genes = adata.var_names.astype(str).to_numpy()
    labels = adata.obs['perturbation'].astype(str).to_numpy()
    ctrl_mask = adata.obs['is_control'].astype(bool).to_numpy()

    granularity = args.granularity
    if granularity == 'donor_env' and 'donor_id' not in adata.obs:
        print('[shoko-v2] no donor_id in obs -> fallback to env granularity', flush=True)
        granularity = 'env'

    self_env_rows = {}
    if granularity == 'donor_env':
        donors = adata.obs['donor_id'].astype(str).to_numpy()
        weights = (adata.obs['cell_count'].to_numpy(dtype=np.float64)
                   if 'cell_count' in adata.obs else np.ones(x.shape[0]))
        r, envs, n_ctrl, row_donors = build_donor_env_response_fast(
            x, labels, donors, ctrl_mask, weights)
        print(f'[shoko-v2] donor_env: {len(envs)} points '
              f'({len(set(row_donors))} donors), {n_ctrl} control wells', flush=True)
    else:
        r, envs, n_ctrl = build_env_response(x, labels, ctrl_mask)
        print(f'[shoko-v2] env: {len(envs)} environments, {n_ctrl} control wells', flush=True)

    floor, dup_excess = find_exclusions(x, genes, args.floor_std)
    excluded = floor | dup_excess
    keep = np.array([g not in excluded for g in genes])
    genes_k = genes[keep]

    tf_all = pd.read_csv(args.tf_all, header=None)[0].astype(str).to_numpy()
    tf_set = set(tf_all)
    tfs_present = [g for g in genes_k if g in tf_set]
    tf_cols = np.flatnonzero(np.isin(genes_k, tfs_present))
    print(f'[shoko-v2] {len(tfs_present)} TFs present, {len(excluded)} genes excluded', flush=True)

    y = r[:, keep]
    y = y - y.mean(axis=0)

    env_set = set(envs)
    env_row = {}
    for i, e in enumerate(envs):
        env_row.setdefault(e, []).append(i)  # donor_env: multiple rows per env
    self_env_rows = {i: env_row[t] for i, t in enumerate(tfs_present) if t in env_set}
    if args.signal == 'ridge' and granularity == 'donor_env' and self_env_rows:
        # ridge + donor_env + self-perturbed TFs: row-downdate not implemented
        raise NotImplementedError('ridge + donor_env with self-perturbed TFs')

    z = np.asarray(x[:, keep], dtype=np.float64)
    z = (z - z.mean(axis=0)) / z.std(axis=0)

    if args.signal == 'pairwise':
        sm = score_matrix(y, tf_cols, self_env_rows)
        alpha_sel = None
    else:
        sm, alpha_sel = ridge_scores(y, tf_cols)
    for i in range(len(tfs_present)):
        sm[i, tf_cols[i]] = -np.inf

    if args.fuse_expr == 'rank50':
        em = expr_corr_matrix(z, tf_cols)
        for i in range(len(tfs_present)):
            em[i, tf_cols[i]] = -np.inf
        sm = rank_fuse(sm, em)
        print('[shoko-v2] rank50 fusion applied', flush=True)

    cand = extract_topn(sm, tfs_present, tf_cols, genes_k, args.top_n_candidates,
                        z, args.pair_corr_max, args.pool_factor)
    print(f'[shoko-v2] candidate ranking: {len(cand)} edges', flush=True)

    total = int(args.max_n_links)
    base = cand.head(total)
    net, info = build_variant(base, cand, int(args.quota), int(args.layer_cap),
                              total, layer_weight=args.layer_weight)
    print(f"[shoko-v2] L_sat={info['L_sat']} -> L={info['L']}, core={info['core_edges']}, "
          f"sources={info['sources']}", flush=True)

    if int(args.tf_out_cap) > 0:
        k = int(args.tf_out_cap)
        before = len(net)
        net = (net.sort_values('weight', ascending=False, kind='stable')
                  .groupby('source', sort=True).head(k)
                  .sort_values('weight', ascending=False, kind='stable')
                  .reset_index(drop=True))
        info['tf_out_cap'] = k
        info['edges_after_cap'] = int(len(net))
        print(f'[shoko-v2] tf_out_cap={k}: {before} -> {len(net)} edges', flush=True)

    cand.to_csv(args.out_prefix + '_candidates.csv', index=False)
    net.to_csv(args.out_prefix + '_ranking.csv', index=False)

    net_str = net.astype(str)
    out = ad.AnnData(
        X=None,
        uns={
            'method_id': 'shoko',
            'dataset_id': adata.uns['dataset_id'],
            'prediction': net_str[['source', 'target', 'weight']],
        },
    )
    out.write_h5ad(args.prediction)

    summary = {
        'dataset': args.dataset,
        'granularity': granularity,
        'signal': args.signal,
        'layer_weight': args.layer_weight,
        'fuse_expr': args.fuse_expr,
        'n_points': len(envs),
        'n_control_wells': n_ctrl,
        'n_tfs_present': len(tfs_present),
        'n_genes_excluded': int(len(excluded)),
        'candidate_edges': int(len(cand)),
        **info,
        'alpha_summary': (None if alpha_sel is None else {
            'min': float(alpha_sel.min()), 'median': float(np.median(alpha_sel)),
            'max': float(alpha_sel.max())}),
        'runtime_s': round(time.time() - t0, 1),
    }
    with open(args.out_prefix + '_summary.json', 'w') as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1), flush=True)


if __name__ == '__main__':
    main()
