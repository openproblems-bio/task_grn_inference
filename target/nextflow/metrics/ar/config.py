"""
This file centralizes the grouping specifications used across metrics.
"""

METHODS = ['positive_control', 'pearson_corr', 'grnboost', 'ppcor', 'portia', 'scenic', 'geneformer', 'scgpt', 'ppcor', 'scenicplus', 'celloracle', 'figr', 'granie', 'scglue', 'scprint',  'negative_control']

DATASET_GROUPS = {
    "op": {
        "match": ["plate_name", "donor_id", "cell_type", "well"],
        "loose_match": ["plate_name", "donor_id", "cell_type"],
        'anchors': ['donor_id', 'plate_name'], # Anchors are preanalytical/technical variables (cell type should not be an anchor)
        "cv": ["perturbation", "cell_type"],
        "rc_tf_ac": ["perturbation", "cell_type"]
    },
    "parsebioscience": {
        "match": ["donor_id", "cell_type", "well"],
        "loose_match": ["donor_id", "cell_type"],
        'anchors': ["donor_id"],
        "cv": ["perturbation", "cell_type"],
        "rc_tf_ac": ["perturbation", "cell_type"]
    },
    "300BCG": {
        "match": ["donor_id", "cell_type"],
        "loose_match": ["cell_type"],  
        'anchors': ["donor_id"],
        "cv": ["perturbation", "cell_type"],
        "rc_tf_ac": ["perturbation", "cell_type"]
    },
    "ibd_UC": {
        'anchors': ['donor_id'],
        "match": ["donor_id", "cell_type"],
        "loose_match": ["donor_id", "cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
    "ibd_CD": {
        'anchors': ['donor_id'],
        "match": ["donor_id", "cell_type"],
        "loose_match": ["donor_id", "cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
    "replogle": {
        "match": ["perturbation"],
        "loose_match": ["perturbation"],
        "cv": ["perturbation"],
    },
    "xaira_HEK293T": {
        "match": ["perturbation"],
        "loose_match": ["perturbation"],
        "cv": ["perturbation"],
    },
    "xaira_HCT116": {
        "match": ["perturbation"],
        "loose_match": ["perturbation"],
        "cv": ["perturbation"],
    },
    "nakatake": {
        "match": ["perturbation"],
        "loose_match": ["perturbation"],
        "cv": ["perturbation"],
    },
     "norman": {
        "match": ["perturbation"],
        "loose_match": ["perturbation"],
        "cv": ["perturbation"],
    },
     "adamson": {
        "match": ["perturbation"],
        "loose_match": ["perturbation"],
        "cv": ["perturbation"],
    },

}

DATASETS = list(DATASET_GROUPS.keys())

DATASETS_CELLTYPES = {
    "replogle": "K562",
    "norman": "K562",
    "adamson": "K562",
    "xaira_HEK293T": "HEK293T",
    "xaira_HCT116": "HCT116",
    "op": "PBMC",
    "parsebioscience": "PBMC",
    "300BCG": "PBMC",
    "ibd_UC": "PBMC",
    "ibd_CD": "PBMC",
    "nakatake": ""
}

DATASETS_METRICS = {
    'replogle': ['regression', 'ws_distance', 'tf_recovery', 'tf_binding', 'sem', 'gs_recovery'],
    'adamson': ['regression', 'ws_distance', 'tf_binding', 'sem', 'gs_recovery'],
    'norman': ['regression', 'ws_distance', 'tf_binding', 'sem', 'gs_recovery'],
    'nakatake': ['regression', 'sem', 'gs_recovery'],
    'op': ['regression', 'vc', 'rc_tf_act', 'tf_binding', 'sem', 'rc', 'anchor_regression', 'gs_recovery'],
    '300BCG': ['regression', 'vc', 'rc_tf_act', 'tf_binding', 'sem', 'rc', 'anchor_regression', 'gs_recovery'],
    'ibd_UC': ['regression', 'vc', 'tf_binding', 'sem', 'anchor_regression', 'gs_recovery'],
    'ibd_CD': ['regression', 'vc', 'tf_binding', 'sem', 'anchor_regression', 'gs_recovery'],
    'parsebioscience': ['regression', 'vc', 'rc_tf_act', 'tf_binding', 'sem', 'rc', 'anchor_regression', 'gs_recovery'],
    'xaira_HEK293T': ['regression', 'ws_distance', 'tf_recovery', 'tf_binding', 'sem', 'gs_recovery'],
    'xaira_HCT116': ['regression', 'ws_distance', 'tf_recovery', 'tf_binding', 'sem', 'gs_recovery'],
}


METRICS_DATASETS = {}
for dataset, metrics in DATASETS_METRICS.items():
    for metric in metrics:
        METRICS_DATASETS.setdefault(metric, []).append(dataset)

ORDERED_METRICS = [
       'r2_precision', 'r2_recall', 
       'ws_precision', 'ws_recall', 
       'vc', 
       'sem', 
       't_rec_precision', 't_rec_recall', 
       'tfb_precision', 'tfb_recall', 
       'rc_tf_act',
       'anchor_regression'
       ]
# ORDERED_METRICS_C = [c for c in ORDERED_METRICS]

def generate_config_env(output_path='src/utils/dataset_config.env'):
    """Generate a simple env-style config file with dataset-specific configurations."""
    
    with open(output_path, 'w') as f:
        f.write("# Auto-generated dataset configuration\n")
        f.write("# Format: DATASET_VARIABLE=value\n\n")
        
        # Cell types
        f.write("# Cell types\n")
        for dataset, cell_type in DATASETS_CELLTYPES.items():
            var_name = f"CELLTYPE_{dataset}"
            f.write(f'{var_name}="{cell_type}"\n')
        
        # Metrics
        f.write("\n# Metrics (comma-separated)\n")
        for dataset, metrics in DATASETS_METRICS.items():
            var_name = f"METRICS_{dataset}"
            metrics_str = ",".join(metrics)
            f.write(f'{var_name}="{metrics_str}"\n')
    
    print(f"Config file generated at: {output_path}")
    return output_path


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate dataset configuration')
    parser.add_argument('--output', type=str, default='src/utils/dataset_config.env',
                       help='Output path for the config file')
    args = parser.parse_args()
    
    generate_config_env(args.output)

