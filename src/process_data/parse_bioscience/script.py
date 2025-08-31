import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os

import pandas as pd
import os
import argparse
import sys

argument_parser = argparse.ArgumentParser(description='Process Replogle data for GRN inference.')
argument_parser.add_argument('--run_test', action='store_true', help='Run in test mode with a subset of data.')
args = argument_parser.parse_args() 



## VIASH START
par =  {'run_test': args.run_test}
## VIASH END

meta = {
    'resources_dir': 'src/process_data/'
}
sys.path.append(meta["resources_dir"])

from helper_data import wrapper_large_perturbation_data


def split_data_perturbation(adata: ad.AnnData, train_share=.5):
    obs = adata.obs
    
    unique_perts = obs['perturbation'].unique()

    control_perturbs = obs[obs['is_control']]['perturbation'].unique()
    n_train_perts = int(len(unique_perts) * train_share)

    # sample for train
    np.random.seed(32)
    train_perturbs = np.random.choice(unique_perts, size=n_train_perts, replace=False)
    test_perturbs = np.setdiff1d(unique_perts, train_perturbs)

    train_perturbs = np.concatenate([train_perturbs, control_perturbs])
    test_perturbs = np.concatenate([test_perturbs, control_perturbs])

    print(f"Train total: {len(train_perturbs)}, Test total: {len(test_perturbs)}")

    return train_perturbs.astype(str), test_perturbs.astype(str)
    
def add_metadata(adata):
    adata.uns['dataset_summary'] = '10 Million Human PBMCs in a Single Experiment'
    adata.uns['dataset_description'] = """Cryopreserved PBMCs from twelve healthy donors were purchased from a commercial vendor. Samples were thawed in a 37༠C water bath, transferred to a 50 mL centrifuge tube, diluted dropwise with warm FBS media, centrifuged, and washed with cold FBS. All samples had a viability >90% after thawing.

                    For each donor, cells were seeded at 1 million cells per well in a 96-well plate for a total of 12 plates across all donors. Cells were treated with 90 different cytokines or PBS (control) for 24 hours resulting in a total of 1,092 experimental conditions. Next, cells were transferred into a deep-well plate, washed with PBS and centrifuged. Once supernatant was removed, cells were fixed with the Evercode Cell Fixation v3 High Throughput Plate-Based Workflow adopted for use with an Integra Assist Plus instrument. In the end, fixed samples were aliquoted into PCR plates, and stored at -80༠C. On a day prior to running the downstream Evercode whole transcriptome experiment, aliquots of fixed samples were thawed in a 37༠C water bath and counted in batches.

                    Fixed samples were processed with Evercode WT Mega v3. After barcoding, 62.45% of cells were retained. Sequencing libraries were made using a Hamilton liquid handler instrument. Libraries were sequenced using the Ultima Genomics sequencing platform and achieved ~31,000 mean reads per cell.

                    After demultiplexing, sequencing data were processed with the Parse Analysis Pipeline v1.4.0. Data were integrated with Harmony, filtered to remove low quality cells, cell types classified with Azimuth, and annotations finalized manually.'
                    """
    adata.uns['data_reference'] = "10 Million Human PBMCs in a Single Experiment, https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/; Parse Biosciences; Accessed 2025 August 8"
    adata.uns['data_url'] = 'https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/'
    adata.uns['dataset_id'] = 'parsebioscience'
    adata.uns['dataset_name'] = 'Parse bioscience (PBMC cytokines)'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    return adata
def format_adata(adata):
    # format
    adata.obs['is_control'] = adata.obs['treatment'] == 'PBS'
    adata.obs = adata.obs[['group', 'cell_type', 'cytokine', 'donor', 'is_control']]
    adata.obs = adata.obs.rename({'donor': 'donor_id'}, axis=1)
    adata.obs['cell_type_minor'] = adata.obs['cell_type']
    cell_type_map = {
        'B Intermediate/Memory': 'B',
        'B Naive': 'B',
        'CD14 Mono': 'MONO',
        'CD16 Mono': 'MONO',
        'CD4 Memory': 'CD4T',
        'CD4 Naive': 'CD4T',
        'Treg': 'CD4T',
        'CD8 Memory': 'CD8T',
        'CD8 Naive': 'CD8T',
        'NK': 'NK',
        'NK CD56bright': 'NK',
        'NKT': 'NK',
        'MAIT': 'CD8T',  # can also be CD4/CD8 double negative, but often grouped under CD8T
        'ILC': 'NK',     # grouped under innate lymphoid/NK for many atlases
        'Plasmablast': 'B',
        'HSPC': None,    # does not map cleanly to any of the 5 requested groups
        'cDC': None,     # dendritic cell, not part of B, MONO, CD8T, CD4T, NK
        'pDC': None      # plasmacytoid DC, same
    }
    adata.obs['cell_type'] = adata.obs['cell_type_minor'].map(cell_type_map)
    adata = adata[~adata.obs['cell_type'].isna(), :]
    adata.obs['well'] = adata.obs['group'].str.split('_').str[-1]
    adata.obs.rename(columns={'cytokine': 'perturbation'}, inplace=True)
    adata.obs['perturbation_type'] = 'cytokine'
    return adata
def main():
    print('Loading data of parsebioscience', flush=True)
    adata = ad.read_h5ad('/vol/projects/CIIM/perturbation_data/Parse_10M_PBMC_cytokines.h5ad', backed='r')
    group_keys = ['cell_type', 'cytokine', 'donor', 'bc1_well']
    for key in group_keys:
        adata.obs[key] = adata.obs[key].astype('str')

    adata.obs['group'] = adata.obs[group_keys].astype(str).agg('_'.join, axis=1)

    if par['run_test']:  # to test
       
        print('Using test data', flush=True)
        # Select one cell per group
        cell_indices = adata.obs.groupby('group').apply(lambda x: x.index[0]).values

        # Create a new AnnData object in memory (small)
        adata_subset = adata[cell_indices, :].to_memory()    
        adata = adata_subset
        cell_count_t = 1
        
    else:
        print('Using full data', flush=True)
        adata = adata.to_memory()
        cell_count_t = 10
    # QC
    min_genes = 10
    min_cell = adata.obs['bc1_well'].nunique()*10

    adata = adata[(adata.obs['gene_count']>min_genes) & (adata.obs['gene_count']<5000), adata.var['n_cells']>min_cell]

    # format the adata
    adata = format_adata(adata)


    # send to main function
    wrapper_large_perturbation_data(adata, split_func=split_data_perturbation,
        covariates=['cell_type_minor', 'perturbation', 'donor_id', 'well'], 
        add_metadata=add_metadata,
        save_name='parsebioscience',
        qc_perturbation_effect=False
        )


if __name__ == '__main__':
    main()


    