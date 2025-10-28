datasets_metrics = {
    'replogle': ['regression', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'adamson': ['regression', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'norman': ['regression', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'nakatake': ['regression', 'sem'],
    'op': ['regression',  'sem',  'tf_binding', 'replica_consistency'],
    '300BCG': ['regression', 'sem',  'tf_binding', 'replica_consistency'],
    'ibd': ['regression', 'sem',  'tf_binding', 'replica_consistency'],
    'parsebioscience': ['regression','sem', 'tf_binding', 'replica_consistency'],
    'xaira_HEK293T': ['regression', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'xaira_HCT116': ['regression', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
}

metrics_datasets = {}
for dataset, metrics in datasets_metrics.items():
    for metric in metrics:
        metrics_datasets.setdefault(metric, []).append(dataset)

