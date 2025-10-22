from typing import List

import numpy as np


class WeightedDegree(object):

    def __init__(self, node_idx: int, degree: int, weight: float):
        """Class used to sort nodes by in-degree or out-degree.

        Ties are broken by using the weights from the GRN.

        Args:
            node_idx
            degree: in-degree or out-degree of a node.
            weight: sum of weights of incoming or outcoming edges.
        """
        self.node_idx: int = int(node_idx)
        self.degree: int = int(degree)
        self.weight: float = float(weight)

    def __lt__(self, other: "WeightedDegree") -> bool:
        if self.degree < other.degree:
            return True
        elif self.degree == other.degree:
            return self.weight < other.weight
        else:
            return False

    def __eq__(self, other: "WeightedDegree") -> bool:
        return (self.degree == other.degree) and (self.weight == other.weight)

    @staticmethod
    def from_grn(A: np.ndarray, incoming: bool = True) -> List["WeightedDegree"]:
        if not incoming:
            A = A.T
        D = (A != 0)
        degrees = np.sum(D, axis=0)
        weights = np.sum(A, axis=0)
        return [WeightedDegree(i, degrees[i], weights[i]) for i in range(len(degrees))]


def create_grn_baseline(A):
    """
    Deterministic baseline for a directed simple graph.
    Preserves in/out degree sequences of A (A may be weighted/signed; nonzeros indicate edges).
    Returns a weighted B: same topology as the deterministic baseline, with the exact weights from A
    reassigned deterministically to different edges (no randomness).
    """

    n = A.shape[0]

    # Order nodes by degrees, with explicit tie-breaking by weight
    in_degrees = WeightedDegree.from_grn(A, incoming=True)
    out_degrees = WeightedDegree.from_grn(A, incoming=False)
    in_order = np.argsort(in_degrees)[::-1]
    out_order = np.argsort(out_degrees)[::-1]

    # Baseline GRN
    B = np.zeros_like(A)

    for u in out_order:
        k = out_degrees[u].degree
        if k == 0:
            continue
        # Greedily fill from current in_order
        picks = []
        for v in in_order:
            if v == u or B[u, v] == 1 or in_degrees[v].degree == 0:
                continue
            picks.append(v)
            if len(picks) == k:
                break

        # Deterministic repair if not enough picks
        if len(picks) < k:
            # Try swaps: reassign some of u's future edges by freeing capacity deterministically
            # Here, we do a deterministic second pass allowing one-time edge relocation
            for v in in_order:
                if v == u or in_degrees[v].degree == 0 or v in picks:
                    continue

                # Find w already chosen where we can swap capacity
                swapped = False
                for w in picks:
                    # Check if there exists x != u with B[x, w]==1 and B[x, v]==0 to swap (x,w)->(x,v)
                    # Deterministic scan by x id
                    for x in range(n):
                        if x == u: 
                            continue
                        if B[x, w] == 1 and B[x, v] == 0 and x != v and x != w:
                            B[x, w] = 0
                            B[x, v] = 1
                            in_degrees[w].degree += 1
                            in_degrees[v].degree -= 1
                            picks.append(v)
                            swapped = True
                            break
                    if swapped or len(picks) == k:
                        break
                if len(picks) == k:
                    break

        if len(picks) < k:
            raise ValueError("Directed degree sequence not realizable with simple digraph under constraints.")

        # Place edges for u
        for v in picks:
            B[u, v] = 1
            in_degrees[v].degree -= 1
        #out_degrees[v].degree = 0
        out_degrees[u].degree = 0

        # Stable re-sort by residual in-degree, with explicit tie-breaking by weight
        in_order = np.argsort(in_degrees)[::-1]

    # Zero diagonal guarantee
    np.fill_diagonal(B, 0)

    # Recompute initial incoming stats from A
    init_in_degrees = WeightedDegree.from_grn(A, incoming=True)

    # Convenience arrays for deterministic target ranking
    # (higher degree first, then higher total incoming weight, then smaller id)
    init_in_deg_arr = np.array([wd.degree for wd in init_in_degrees])
    init_in_wgt_arr = np.array([wd.weight for wd in init_in_degrees])

    for u in range(n):

        # Outgoing weights in A
        mask_A_u = (A[u, :] != 0)
        if not np.any(mask_A_u):
            continue
        orig_targets = np.where(mask_A_u)[0]
        W = A[u, orig_targets].astype(float)

        # Sort weights deterministically: by |w| desc, then w asc, then orig target id asc
        # (lexsort uses last key as primary)
        w_keys_3 = np.abs(W) * -1.0       # primary: |w| descending
        w_keys_2 = W                      # secondary: actual value ascending (keeps sign order stable)
        w_keys_1 = orig_targets           # tertiary: original target id ascending
        order_w = np.lexsort((w_keys_1, w_keys_2, w_keys_3))
        W_sorted = W[order_w]

        # Targets in B for this source, ranked by original A's incoming difficulty/salience
        # Rank by: in-degree desc, then incoming-weight desc, then id asc
        targets_B = np.where(B[u, :] == 1)[0]
        if targets_B.size == 0:
            continue
        t_deg = init_in_deg_arr[targets_B]
        t_wgt = init_in_wgt_arr[targets_B]
        t_id  = targets_B
        order_t = np.lexsort((t_id, -t_wgt, -t_deg))  # primary: -t_deg, then -t_wgt, then t_id
        targets_ranked = targets_B[order_t]

        # Sanity: degrees should match. If not, deterministically trim/pad.
        kA = W_sorted.size
        kB = targets_ranked.size
        if kA > kB:
            # Drop evenly from both ends to avoid bias
            excess = kA - kB
            left = excess // 2
            right = excess - left
            W_sorted = W_sorted[left: kA - right]
        elif kA < kB:
            W_sorted = np.concatenate([W_sorted, np.zeros(kB - kA, dtype=float)], axis=0)

        # Assign exact weights to new edges deterministically
        B[u, targets_ranked] = W_sorted

    return B
