"""
This file centralizes the grouping specifications used across metrics.
"""

DATASET_GROUPS = {
    "op": {
        "match": ["plate_name", "donor_id", "cell_type", "well"],
        "loose_match": ["plate_name", "donor_id", "cell_type"],
        'anchors': ['donor_id', 'plate_name'], # Anchors are preanalytical/technical variables (cell type should not be an anchor)
        "cv": ["perturbation", "cell_type"],
    },
    "parsebioscience": {
        "match": ["donor_id", "cell_type", "well"],
        "loose_match": ["donor_id", "cell_type"],
        'anchors': ["donor_id"],
        "cv": ["perturbation", "cell_type"],
    },
    "300BCG": {
        "match": ["donor_id", "cell_type"],
        "loose_match": ["cell_type"],  
        'anchors': ["donor_id"],
        "cv": ["perturbation", "cell_type"],
    },
    "IBD": {
        'anchors': ['donor_id'],
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

