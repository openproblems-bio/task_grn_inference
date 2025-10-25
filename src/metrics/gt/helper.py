from typing import Optional, Dict

import numpy as np
import pandas as pd
from sklearn.metrics import recall_score, precision_score, matthews_corrcoef, roc_auc_score, average_precision_score
from util import read_prediction


# For reproducibility
seed = 42
np.random.seed(seed)


def load_grn_from_dataframe(
        df: pd.DataFrame,
        tf_dict: Dict[str, int],
        tg_dict: Dict[str, int]
) -> np.ndarray:
    G = np.zeros((len(tf_dict), len(tg_dict)), dtype=float)
    for tf_name, tg_name, weight in zip(df['source'], df['target'], df['weight']):
        if tf_name not in tf_dict:
            continue
        if tg_name not in tg_dict:
            continue
        i = tf_dict[tf_name]
        j = tg_dict[tg_name]
        G[i, j] = weight
    return G


def main(par):

    # Load ground-truth edges (consider only TFs listed in the file loaded hereabove)
    true_graph = pd.read_csv(par['ground_truth'])
    if not {'Gene1', 'Gene2'}.issubset(set(true_graph.columns)):
        raise ValueError("ground_truth must have columns: 'Gene1', 'Gene2'")
    if 'weight' not in true_graph:
        true_graph['weight'] = np.ones(len(true_graph))
    true_graph.rename(columns={"Gene1": "source", "Gene2": "target"}, copy=False, inplace=True)

    # Load inferred GRN
    prediction = read_prediction(par)
    assert prediction.shape[0] > 0, 'No links found in the network'
    if not {'source', 'target', 'weight'}.issubset(set(prediction.columns)):
        raise ValueError("prediction must have columns: 'source', 'target', 'weight'")

    # Intersect TF lists, intersect TG lists
    tf_names = set(true_graph['source'].unique())  #.intersection(set(tf_all))
    tg_names = set(true_graph['target'].unique())
    #tf_names = set(true_graph['source'].unique()).intersection(set(prediction['source'].unique()))  #.intersection(set(tf_all))
    #tg_names = set(true_graph['target'].unique()).intersection(set(prediction['target'].unique()))
    tf_dict = {tf_name: i for i, tf_name in enumerate(tf_names)}
    tg_dict = {tg_name: i for i, tg_name in enumerate(tg_names)}

    # Reformat GRNs as NumPy arrays
    A = load_grn_from_dataframe(true_graph, tf_dict, tg_dict)
    G = load_grn_from_dataframe(prediction, tf_dict, tg_dict)
    G = np.abs(G)

    # Evaluate inferred GRN
    tf_binding_precision = precision_score((A != 0).flatten(), (G != 0).flatten())
    tf_binding_recall = recall_score((A != 0).flatten(), (G != 0).flatten())
    final_score = roc_auc_score(A.flatten(), G.flatten())

    summary_df = pd.DataFrame([{
        'tf_binding_precision': tf_binding_precision,
        'tf_binding_recall': tf_binding_recall,
        'final_score': final_score
    }])

    return summary_df
