"""
This file centralizes the grouping specifications used across metrics.
"""

DATASET_GROUPS = {
    "op": {
        "match": ["plate_name", "donor_id", "cell_type", "well"],
        "loose_match": ["donor_id", "cell_type"],
        'anchors': ['donor_id', 'plate_name', 'cell_type'], #should cell type be an anchor? 
        "cv": ["perturbation", "cell_type"],
    },
    "parsebioscience": {
        "match": ["donor_id", "cell_type", "well"],
        "loose_match": ["donor_id", "cell_type"],
        'anchor': ["donor_id", "cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
    "300BCG": {
        "match": ["donor_id", "cell_type"],
        "loose_match": ["cell_type"],
        'anchor': ["donor_id", "cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
}