#!/usr/bin/env python3
"""eval23 stream_core.py — streaming sufficient statistics for the final GRN
method (依頼23, F2/F3 of prereg_simplify.md).

Computes, in ONE chunked pass over cells (never materializing dense n x G):
  - per-gene sum / sumsq          -> mean, std, exclusions
  - per-group (env or donor_env) response sums + counts -> response profiles
  - TF x G cross-product (float64 accumulator) -> exact cell-level Pearson (em)
Optional streaming randomized PCA (for the W retention gate) without
holding z: chunked sketch + QR-normalized power iterations + chunked
B = Q^T Z. Same algorithm family as sklearn randomized_svd (approximate
top-k PCs; used only for the retention gate / W reselection).

All statistics are mathematically identical to the current dense pipeline
(F2) and the pool filter becomes an em lookup (F3).
"""
import h5py
import numpy as np
import scipy.sparse as sp


class H5AdReader:
    """Minimal streaming reader for h5ad inference data.

    anndata read_h5ad(backed='r') eagerly loads LAYERS into RAM (only X is
    truly backed) — 45GB peak on the 1M-cell synthetic file. This reader
    pulls obs/var metadata and layer row-chunks directly via h5py, so peak
    memory is O(chunk), independent of n_cells.
    """

    def __init__(self, path, layer_key):
        self.f = h5py.File(path, 'r')
        var = self.f['var']
        idx_key = var.attrs.get('_index', '_index')
        raw = var[idx_key][:]
        self.genes = np.array([g.decode() if isinstance(g, bytes) else str(g)
                               for g in raw])
        node = self.f[f'layers/{layer_key}']
        if isinstance(node, h5py.Group):  # csr_matrix
            self.kind = 'csr'
            self.data = node['data']
            self.indices = node['indices']
            self.indptr = node['indptr']
            self.shape = tuple(int(v) for v in node.attrs['shape'])
        else:  # dense
            self.kind = 'dense'
            self.dset = node
            self.shape = (int(node.shape[0]), int(node.shape[1]))

    def has_obs(self, key):
        return key in self.f['obs']

    def obs_col(self, key):
        node = self.f['obs'][key]
        if isinstance(node, h5py.Group):  # categorical
            cats = node['categories'][:]
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                             for c in cats])
            codes = node['codes'][:]
            if (codes < 0).any():
                raise ValueError(f'obs/{key}: NaN categorical codes not supported')
            return cats[codes]
        return node[:]

    def layer_rows(self, s, e):
        """Dense float32 (e-s, G) chunk read straight from disk."""
        if self.kind == 'dense':
            return np.asarray(self.dset[s:e], dtype=np.float32)
        a = int(self.indptr[s])
        b = int(self.indptr[e])
        dat = self.data[a:b]
        idx = self.indices[a:b].astype(np.int64, copy=False)
        ip = self.indptr[s:e + 1].astype(np.int64) - a
        m = sp.csr_matrix((dat, idx, ip), shape=(e - s, self.shape[1]))
        return m.toarray().astype(np.float32, copy=False)

    def close(self):
        self.f.close()


class StreamStats:
    """Chunked accumulation over cells.

    Parameters
    ----------
    genes_k : array of gene names (post-exclusion is applied later; here full)
    tf_cols : int array, columns of TFs in the gene axis
    group_keys : array-like of per-cell group labels (env or donor_env strings)
    ctrl_mask : bool array per cell
    chunk : cells per chunk
    """

    def __init__(self, n_genes, tf_cols, group_keys, ctrl_mask, chunk=25_000,
                 donors=None, weights=None):
        self.G = n_genes
        self.tf_cols = np.asarray(tf_cols)
        self.chunk = chunk
        self.sum = np.zeros(n_genes, dtype=np.float64)
        self.sumsq = np.zeros(n_genes, dtype=np.float64)
        self.n = 0
        # group accumulators
        self.groups = {}
        self.group_keys = np.asarray(group_keys)
        self.ctrl_mask = np.asarray(ctrl_mask)
        uniq = list(dict.fromkeys(self.group_keys))
        self.g2i = {g: i for i, g in enumerate(uniq)}
        self.n_groups = len(uniq)
        self.gsum = np.zeros((self.n_groups, n_genes), dtype=np.float64)
        self.gcount = np.zeros(self.n_groups, dtype=np.float64)
        self.ctrl_sum = np.zeros(n_genes, dtype=np.float64)
        self.n_ctrl = 0
        # optional cell weights (used ONLY for the response matrix, matching
        # build_donor_env_response_fast's cell_count weighting)
        self.weights = np.asarray(weights, dtype=np.float64) if weights is not None else None
        if self.weights is not None:
            self.gsum_w = np.zeros((self.n_groups, n_genes), dtype=np.float64)
            self.gw = np.zeros(self.n_groups, dtype=np.float64)
        # optional per-donor control (donor_env granularity: pb path)
        self.donors = np.asarray(donors) if donors is not None else None
        if self.donors is not None:
            self.d2i = {d: i for i, d in enumerate(dict.fromkeys(self.donors))}
            self.dctrl_sum = np.zeros((len(self.d2i), n_genes), dtype=np.float64)
            self.dctrl_n = np.zeros(len(self.d2i), dtype=np.float64)
            # donor of each group key "donor|env"
            self.group_donor = np.array([self.d2i[g.split('|')[0]] for g in uniq])
        # TF x G cross-product (raw, uncentered, unweighted)
        self.xtx_tf = np.zeros((len(self.tf_cols), n_genes), dtype=np.float64)

    def update(self, x_chunk, rows):
        """x_chunk: (m, G) float32 dense chunk (CSR chunks densified per chunk
        only). rows: global row indices of the chunk."""
        m = x_chunk.shape[0]
        xd = x_chunk.astype(np.float64)
        self.sum += xd.sum(axis=0)
        self.sumsq += (xd ** 2).sum(axis=0)
        self.n += m
        gidx = np.array([self.g2i[g] for g in self.group_keys[rows]])
        for gi in np.unique(gidx):
            sel = gidx == gi
            self.gsum[gi] += xd[sel].sum(axis=0)
            self.gcount[gi] += sel.sum()
            if self.weights is not None:
                w = self.weights[rows][sel][:, None]
                self.gsum_w[gi] += (xd[sel] * w).sum(axis=0)
                self.gw[gi] += w.sum()
        cm = self.ctrl_mask[rows]
        if cm.any():
            self.ctrl_sum += xd[cm].sum(axis=0)
            self.n_ctrl += int(cm.sum())
            if self.donors is not None:
                for d in np.unique(self.donors[rows][cm]):
                    sel = cm & (self.donors[rows] == d)
                    if self.weights is not None:
                        w = self.weights[rows][sel][:, None]
                        self.dctrl_sum[self.d2i[d]] += (xd[sel] * w).sum(axis=0)
                        self.dctrl_n[self.d2i[d]] += w.sum()
                    else:
                        self.dctrl_sum[self.d2i[d]] += xd[sel].sum(axis=0)
                        self.dctrl_n[self.d2i[d]] += sel.sum()
        self.xtx_tf += xd[:, self.tf_cols].T @ xd

    # --------------------------------------------------------- derived ----
    def mean_std(self):
        mean = self.sum / self.n
        var = self.sumsq / self.n - mean ** 2
        var[var < 0] = 0.0
        return mean, np.sqrt(var)

    def response_matrix(self):
        """env response = group mean - control mean (matches
        build_env_response: unweighted means, global control).
        With donors (+weights): cell_count-weighted group means minus
        per-donor weighted control mean (matches build_donor_env_response_fast)."""
        if self.donors is None:
            ctrl_mean = self.ctrl_sum / self.n_ctrl
            r = self.gsum / self.gcount[:, None] - ctrl_mean[None, :]
        else:
            dmean = self.dctrl_sum / np.maximum(self.dctrl_n, 1e-12)[:, None]
            if self.weights is not None:
                r = self.gsum_w / np.maximum(self.gw, 1e-12)[:, None] - dmean[self.group_donor]
            else:
                r = self.gsum / self.gcount[:, None] - dmean[self.group_donor]
        return r.astype(np.float64)

    def em_signed(self):
        """Exact signed Pearson TF x G from the cross-product."""
        mean, std = self.mean_std()
        std_safe = np.where(std == 0, 1.0, std)
        cov = self.xtx_tf / self.n - np.outer(mean[self.tf_cols], mean)
        r = cov / np.outer(std_safe[self.tf_cols], std_safe)
        return r


def streaming_randomized_pca(load_chunk, n_cells, n_genes, mean, std, n_components=5,
                             n_oversamples=10, n_iter=3, seed=0):
    """Streaming randomized PCA on implicitly standardized data.

    load_chunk(start, stop) -> (m, G) float32 RAW chunk (unstandardized).
    Standardization applied on the fly with precomputed mean/std.
    Returns (U, S, Vt) like sklearn's randomized_svd (approximate).

    Passes over data: 1 (sketch) + 2*n_iter (power iterations, QR-normalized)
    + 1 (B = Q^T Z). n_iter=3 reproduces sklearn randomized_svd
    (n_iter='auto' -> 7) retention medians to ~3e-3 on pb/nakatake
    (test_pca_power.py); n_iter=0 is NOT sufficient (subspace contaminated
    by the spectral tail).
    """
    rng = np.random.default_rng(seed)
    k = n_components + n_oversamples
    omega = rng.standard_normal((n_genes, k))
    std_safe = np.where(std == 0, 1.0, std)
    CH = 10_000  # rows per chunk: 10k x 20k float64 = 1.6GB transient (100k OOMs at 1M cells)

    def z_matvec(A):  # Y = Z @ A, one chunked pass
        Y = np.zeros((n_cells, A.shape[1]), dtype=np.float64)
        for s in range(0, n_cells, CH):
            e = min(s + CH, n_cells)
            zc = load_chunk(s, e).astype(np.float64)
            zc -= mean
            zc /= std_safe
            Y[s:e] = zc @ A
        return Y

    def zt_matvec(Q):  # A = Z^T @ Q, one chunked pass
        A = np.zeros((n_genes, Q.shape[1]), dtype=np.float64)
        for s in range(0, n_cells, CH):
            e = min(s + CH, n_cells)
            zc = load_chunk(s, e).astype(np.float64)
            zc -= mean
            zc /= std_safe
            A += zc.T @ Q[s:e]
        return A

    Q, _ = np.linalg.qr(z_matvec(omega))
    for _ in range(n_iter):
        Q, _ = np.linalg.qr(z_matvec(zt_matvec(Q)))
    Bt = zt_matvec(Q)  # G x k (= Z^T Q; B = Q^T Z is its transpose)
    Vb, S, Ut = np.linalg.svd(Bt, full_matrices=False)
    U = Q @ Ut.T
    return U[:, :n_components], S[:n_components], Vb.T[:n_components]
