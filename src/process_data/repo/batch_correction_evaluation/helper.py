
import scib
import lightgbm as lgb
import pandas as pd
from sklearn.model_selection import cross_validate
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import r2_score, make_scorer, accuracy_score


def run_scib(bulk_adata, layer='lognorm', layer_baseline='n_counts', batch_key='plate_name', label_key='cell_type'):
    bulk_adata.X = bulk_adata.layers[layer_baseline].copy()

    bulk_adata_c = bulk_adata.copy()
    bulk_adata_c.X = bulk_adata_c.layers[layer].copy()

    scib.pp.reduce_data(
        bulk_adata_c, n_top_genes=None, batch_key=batch_key, pca=True, neighbors=True
    )
    rr = scib.metrics.metrics(bulk_adata, bulk_adata_c, batch_key, label_key, organism='human', 
                            # biological conservation (label)
                            nmi_=True, 
                            ari_=False,
                            silhouette_=True,
                            isolated_labels_f1_=False, # there is no isolated cell type
                            isolated_labels_asw_=False, # there is no isolated cell type
                            # biological conservation (label free)
                            cell_cycle_=True,
                            hvg_score_=False,
                            trajectory_=False,
                            # batch correction
                            pcr_=False, 
                            graph_conn_=False,
                            kBET_=True,
                            ilisi_=False,
                            clisi_=False,
                            # Not sure what they are
                            isolated_labels_=False,  # backwards compatibility
                            n_isolated=None,
                            lisi_graph_=False,

                            verbose = 0
                            )
    rr = rr.dropna().T
    return rr 
def run_classifier(adata, layer, batch_key):
    print('GB classifier')
    model = lgb.LGBMClassifier(silent=True, verbose=-1)
    # model = RidgeClassifier()
    X = adata.layers[layer].copy()
    y = adata.obs[batch_key]
    scoring = {
        'accuracy_score': make_scorer(accuracy_score)
    }
    score = 1 - cross_validate(model, X, y, cv=5, scoring=scoring, return_train_score=False)['test_accuracy_score'].mean()
    
    return pd.DataFrame({'Batch classifier':[score]})
