
import os
import gc

import sys
import yaml
import pickle
import tempfile
import contextlib
import hashlib
import shutil
import requests
import traceback
import zipfile
import subprocess
import gc
import gzip
import tarfile
from urllib.request import urlretrieve
import json
import os



# import flatbuffers
import numpy as np
import scanpy as sc
import pandas as pd
import anndata
import pyranges as pr
import polars
import mudata
import scrublet as scr
from sklearn.manifold import TSNE
import pycisTopic.loom
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.lda_models import evaluate_models, run_cgs_models, run_cgs_models_mallet
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
from pycisTopic.diff_features import impute_accessibility, normalize_scores, find_highly_variable_features, find_diff_features
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.gene_activity import get_gene_activity
from pycisTopic.loom import export_region_accessibility_to_loom, export_gene_activity_to_loom
from pycisTopic.clust_vis import find_clusters, run_umap, run_tsne, plot_metadata, plot_topic, cell_topic_heatmap
from scenicplus.wrappers.run_pycistarget import run_pycistarget
# Bug in pycistopic: Import is missing in pycisTopic.loom,
# so TSNE must be dynamically added to the library's namespace.
setattr(pycisTopic.loom, 'TSNE', TSNE)

os.environ['MALLET_MEMORY'] = '200G'
try:
    sys.path.append(meta['resources_dir'])
except NameError:
    pass

def process_peak(par):
    # Get list of samples (e.g., donors)
    print('Collect list of samples')
    adata_atac = anndata.read_h5ad(par['atac'])
    unique_donor_ids = [s.replace(' ', '_') for s in adata_atac.obs.donor_id.cat.categories]
    print(unique_donor_ids)
    unique_cell_types = [s.replace(' ', '_') for s in adata_atac.obs.cell_type.cat.categories]


    # Create one individual ATAC-seq file per donor
    fragments_dict = {}
    for donor_id in unique_donor_ids:
        filepath = os.path.join(par['atac_dir'], f'{donor_id}.tsv')
        print(f'Create tsv file {filepath}')
        adata_atac.obs.cell_type = [s.replace(' ', '_') for s in adata_atac.obs.cell_type]
        adata_atac.obs.donor_id = [s.replace(' ', '_') for s in adata_atac.obs.donor_id]
        adata_atac_one_donor = adata_atac[adata_atac.obs.donor_id == donor_id]
        coo_matrix = adata_atac_one_donor.X.tocoo()
        rows = coo_matrix.row
        cols = coo_matrix.col
        counts = coo_matrix.data
        row_names = adata_atac_one_donor.obs.index[rows]
        coord_names = adata_atac_one_donor.var.index[cols]
        del  adata_atac_one_donor
        gc.collect()
        d = {
            'chromosome': np.asarray([s.split(':')[0] for s in coord_names]),
            'start': np.asarray([int(s.split(':')[1].split('-')[0]) for s in coord_names], dtype=np.uint32),
            'end': np.asarray([int(s.split(':')[1].split('-')[1]) for s in coord_names], dtype=np.uint32),
            'barcode': [barcode.replace('-', '') for barcode in row_names],
            'count': np.squeeze(np.asarray(counts.astype(np.uint16)))
        }
        df = pd.DataFrame(d, copy=False)
        df.to_csv(filepath, sep='\t', index=False, header=False)

        # Compress tsv file
        compressed_filepath = filepath + '.gz'
        print(f'Sort and compress tsv file {filepath}')

        # Sort file by genomic coordinates
        sorted_filepath = filepath + '.sorted.tsv'
        os.system(f'sort -k1,1 -k2,2n {filepath} > {sorted_filepath}')

        # Decompression
        os.system(f'bgzip -c {sorted_filepath} > {compressed_filepath}')

        fragments_dict[donor_id] = compressed_filepath

        # Index file using tabix
        print(f'Index compressed file {compressed_filepath} using tabix')
        subprocess.run(['tabix', compressed_filepath, '-p', 'bed'])
    
    with open(f"{par['fragments_dict']}", 'w') as f:
        json.dump(fragments_dict, f)
    # Collect cell metadata
    print(f'Collect cell metadata')
    donor_ids = [s.replace(' ', '_') for s in adata_atac.obs.donor_id]
    index = [barcode.replace('-', '') + '-' + donor_id for donor_id, barcode in zip(donor_ids, adata_atac.obs.index)]
    cell_data = pd.DataFrame({
        'cell_type': [s.replace(' ', '_') for s in adata_atac.obs.cell_type.to_numpy()],
        'donor_id': [s.replace(' ', '_') for s in adata_atac.obs.donor_id.to_numpy()]
    }, index=index, copy=False)
    

    # Get chromosome sizes
    # if True:
    #     print(f'Download chromosome sizes from UCSC')
    #     chromsizes = pd.read_table(
    #         'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    #         header=None,
    #         names=['Chromosome', 'End']
    #     )
    # else:
    par['chromsizes'] = par['temp_dir'] + 'chromsizes.tsv'

    response = requests.get('http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes')
    with open(par['chromsizes'], "wb") as file:
        file.write(response.content)

    chromsizes = pd.read_csv(par['chromsizes'], sep='\t', names=['Chromosome', 'End'])

    chromsizes.insert(1, 'Start', 0)

    print(f'Start pseudobulking')
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bed_files'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bw_files'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/MACS'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/tmp'), exist_ok=True)

    # Check if individual fragment files are available for each cell type
    bed_paths = {}
    for cell_type in unique_cell_types:
        filepath = os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bed_files', f'{cell_type}.fragments.tsv.gz')
        bed_paths[cell_type] = filepath


    # Perform pseudobulking for each cell type
    print('Pseudobulk each cell type')
    _, bed_paths = export_pseudobulk(
        input_data=cell_data,
        variable='cell_type',
        sample_id_col='donor_id',
        chromsizes=chromsizes,
        bed_path=os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bed_files'),
        bigwig_path=os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bw_files'),
        path_to_fragments=fragments_dict,
        n_cpu=par['num_workers'],
        temp_dir=os.path.join(par['temp_dir'], 'consensus_peak_calling/tmp'),
        split_pattern='-',
    )
    with open(os.path.join(par['temp_dir'], 'consensus_peak_calling/bed_paths.tsv'), 'wt') as f:
        for v in bed_paths:
            _ = f.write(f'{v}\t{bed_paths[v]}\n')

    # Load paths to pseudobulked samples
    bed_paths = {}
    with open(os.path.join(par['temp_dir'], 'consensus_peak_calling/bed_paths.tsv')) as f:
        for line in f:
            v, p = line.strip().split("\t")
            bed_paths.update({v: p})

    # Call peaks using MACS2
    print('Call peaks using MACS2')
    narrow_peak_dict = peak_calling(
        macs_path = 'macs2',
        bed_paths=bed_paths,
        outdir=os.path.join(os.path.join(par['temp_dir'], 'consensus_peak_calling/MACS')),
        genome_size='hs',
        n_cpu=10,
        input_format='BEDPE',
        shift=73,
        ext_size=146,
        keep_dup='all',
        q_value=0.05
    )

    # Consensus peak calling
    print('Consensus peak calling')
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peak_dict,
        peak_half_width=250,
        chromsizes=chromsizes,
        path_to_blacklist=par['blacklist_path']
    )
    consensus_peaks.to_bed(
        path=os.path.join(par['temp_dir'], 'consensus_peak_calling/consensus_regions.bed'),
        keep=True,
        compression='infer',
        chain=False
    )

    # Download Mallet
    if not os.path.exists(par['MALLET_PATH']):
        url = 'https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz'
        response = requests.get(url)
        with open(os.path.join(par['temp_dir'], 'Mallet-202108-bin.tar.gz'), 'wb') as f:
            f.write(response.content)
        with tarfile.open(os.path.join(par['temp_dir'], 'Mallet-202108-bin.tar.gz'), 'r:gz') as f:
            f.extractall(path=par['temp_dir'])
    del consensus_peaks
    del narrow_peak_dict
    del adata_atac
    gc.collect()
def run_cistopic(par):
    adata_atac = anndata.read_h5ad(par['atac'])
    unique_donor_ids = [s.replace(' ', '_') for s in adata_atac.obs.donor_id.cat.categories]
    print(unique_donor_ids)
    unique_cell_types = [s.replace(' ', '_') for s in adata_atac.obs.cell_type.cat.categories]
    donor_ids = [s.replace(' ', '_') for s in adata_atac.obs.donor_id]
    index = [barcode.replace('-', '') + '-' + donor_id for donor_id, barcode in zip(donor_ids, adata_atac.obs.index)]
    cell_data = pd.DataFrame({
        'cell_type': [s.replace(' ', '_') for s in adata_atac.obs.cell_type.to_numpy()],
        'donor_id': [s.replace(' ', '_') for s in adata_atac.obs.donor_id.to_numpy()]
    }, index=index, copy=False)
    
    with open(f"{par['fragments_dict']}", 'r') as f:
        fragments_dict = json.load(f)
    # run pycistopic
    print('get tss')
    os.makedirs(os.path.join(par['temp_dir'], 'qc'), exist_ok=True)
    subprocess.run([
        'pycistopic', 'tss', 'get_tss',
        '--output', os.path.join(par['temp_dir'], 'qc', 'tss.bed'),
        '--name', 'hsapiens_gene_ensembl',
        '--to-chrom-source', 'ucsc',
        '--ucsc', 'hg38'
    ])

    if par['qc']:  # Whether to perform quality control
        # Compute QC metrics
        print('Perform QC')
        for donor_id in unique_donor_ids:
            filepath = os.path.join(par['atac_dir'], f'{donor_id}.tsv.gz')
            subprocess.run([
                'pycistopic', 'qc',
                '--fragments', os.path.join(par['temp_dir'], 'qc', 'tss.bed'),
                '--regions', os.path.join(par['temp_dir'], 'consensus_peak_calling/consensus_regions.bed'),
                '--tss', os.path.join(par['temp_dir'], 'qc', 'tss.bed'),
                '--output', os.path.join(par['temp_dir'], 'qc', donor_id),
                '--tss_flank_window', '10000',  # Default: 2000
                '--tss_smoothing_rolling_window', '60',  # Default: 10
                '--tss_minimum_signal_window', '5',  # Default: 100
                '--tss_window', '25',  # Default: 50
                '--tss_min_norm', '0.5',  # Default: 0.2
                '--min_fragments_per_cb', '30',  # Default: 10
                '--threads', str(par['num_workers'])
            ])

        # Filter out low quality cells
        print('Filter out low quality cells')
        sample_id_to_barcodes_passing_filters = {}
        sample_id_to_thresholds = {}
        for donor_id in fragments_dict.keys():
            (
                sample_id_to_barcodes_passing_filters[donor_id],
                sample_id_to_thresholds[donor_id]
            ) = get_barcodes_passing_qc_for_sample(
                    sample_id=donor_id,
                    pycistopic_qc_output_dir=os.path.join(par['temp_dir'], 'qc'),
                    unique_fragments_threshold=None, # use automatic thresholding
                    tss_enrichment_threshold=None, # use automatic thresholding
                    frip_threshold=0,
                    use_automatic_thresholds=True,
            )

        # Create cistopic objects
        print('Create cisTopic objects')
        cistopic_obj_list = []
        for donor_id in unique_donor_ids:
            sample_metrics = pl.read_parquet(
                os.path.join(par['temp_dir'], 'qc', donor_id, f'{donor_id}.fragments_stats_per_cb.parquet')
            ).to_pandas().set_index('CB').loc[sample_id_to_barcodes_passing_filters[sample_id]]
            cistopic_obj = create_cistopic_object_from_fragments(
                path_to_fragments=fragments_dict[donor_id],
                path_to_regions=os.path.join(par['temp_dir'], 'consensus_peak_calling/consensus_regions.bed'),
                path_to_blacklist=par['blacklist_path'],
                metrics=sample_metrics,
                valid_bc=sample_id_to_barcodes_passing_filters[sample_id],
                n_cpu = par['num_workers'],
                project=donor_id,
                split_pattern='-'
            )
            cistopic_obj_list.append(cistopic_obj)
    else:

        # Skip QC and create cisTopic objects
        print('Create cisTopic objects without QC')
        cistopic_obj_list = []
        for donor_id in unique_donor_ids:
            cistopic_obj = create_cistopic_object_from_fragments(
                path_to_fragments=fragments_dict[donor_id],
                path_to_regions=os.path.join(par['temp_dir'], 'consensus_peak_calling/consensus_regions.bed'),
                path_to_blacklist=par['blacklist_path'],
                n_cpu=par['num_workers'],
                project=donor_id,
                split_pattern='-'
            )
            cistopic_obj_list.append(cistopic_obj)

    # Add metadata to cistopic objects
    for i in range(len(cistopic_obj_list)):
        cistopic_obj_list[i].add_cell_data(cell_data, split_pattern='-')

    # Infer doublets using scrublet
    # print('Infer doublets using scrublet')
    # for i in range(len(cistopic_obj_list)):
    #     print(f'step 1')
    #     scrub = scr.Scrublet(cistopic_obj_list[i].fragment_matrix.T, expected_doublet_rate=0.1)
    #     print(f'step 2')
    #     doublet_scores, predicted_doublets = scrub.scrub_doublets()
    #     print(f'step 3')
    #     scrub.call_doublets(threshold=0.22)
    #     print(f'step 4')
    #     scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T
    #     print(f'step 5')
    #     cistopic_obj_list[i].add_cell_data(scrublet, split_pattern='-')
    # The code is printing the message "Infer doublets using scrublet completed" in Python.
    # print('Infer doublets using scrublet completed')
    # Combine samples into a single cistopic object
    if len(cistopic_obj_list) == 1:
        cistopic_obj = cistopic_obj_list[0]
    else:
        cistopic_obj = merge(cistopic_obj_list, is_acc=1, split_pattern='-')
    
    # LDA-based topic modeling
    print('Run LDA models', flush=True)
    n_topics = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    # n_topics = [40] #TODO: fix this
    if os.path.exists(par['MALLET_PATH']):
        models = run_cgs_models_mallet(
            cistopic_obj,
            n_topics=n_topics,
            n_cpu=par['num_workers'],
            n_iter=500,
            random_state=555,
            alpha=50,
            alpha_by_topic=True,
            eta=0.1,
            eta_by_topic=False,
            mallet_path=par['MALLET_PATH']
        )
    else:
        print('Could not find Mallet. Running the sequential version of LDA instead.')
        models = run_cgs_models(
            cistopic_obj,
            n_topics=n_topics,
            n_cpu=par['num_workers'],
            n_iter=500,
            random_state=555,
            alpha=50,
            alpha_by_topic=True,
            eta=0.1,
            eta_by_topic=False
        )

    # Model selection
    model = evaluate_models(models, select_model=40, return_model=True)
    cistopic_obj.add_LDA_model(model)

    # Save cistopic objects
    with open(par['cistopic_object'], 'wb') as f:
        pickle.dump(cistopic_obj, f)
        
    del cistopic_obj
    del model
    del cistopic_obj_list
    del adata_atac
    gc.collect()

def process_topics(par):  
    # Load cistopic objects
    with open(par['cistopic_object'], 'rb') as f:
        cistopic_obj = pickle.load(f)
        
    cell_topic = cistopic_obj.selected_model.cell_topic.T
    cell_topic.index = cell_topic.index.str.split('-').str[0]
    cell_topic.to_csv(par['cell_topic'])

    # chromsizes = pd.read_table(os.path.join(par['temp_dir'], "qc", "hg38.chrom_sizes_and_alias.tsv"))
    
    # chromsizes.rename({"# ucsc": "Chromosome", "length": "End"}, axis = 1, inplace = True)
    chromsizes = pd.read_csv(par['chromsizes'], sep='\t', names=['Chromosome', 'End'])
    chromsizes['Start'] = 0
    chromsizes = pr.PyRanges(chromsizes[['Chromosome', 'Start', 'End']])

    # Find clusters
    find_clusters(
        cistopic_obj,
        target='cell',
        k=10,
        res=[0.6, 1.2, 3],
        scale=True,
        split_pattern='-'
    )

    # 2D projections
    run_umap(cistopic_obj, target='cell', scale=True)
    run_tsne(cistopic_obj, target='cell', scale=True)

    # Topic binarization
    region_bin_topics_top_3k = binarize_topics(cistopic_obj, method='ntop', ntop=3_000)
    region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
    binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li')

    # Topic annotation
    topic_annot = topic_annotation(
        cistopic_obj,
        annot_var='cell_type',
        binarized_cell_topic=binarized_cell_topic,
        general_topic_thr=0.2
    )

    # Identify differentially accessible regions
    imputed_acc_obj = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6
    )
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp=0.05,
        min_mean=0.0125,
        max_mean=3,
        max_disp=np.inf,
        n_bins=20,
        n_top_features=None,
        plot=False
    )
    markers_dict = find_diff_features(
        cistopic_obj,
        imputed_acc_obj,
        variable='cell_type',
        var_features=variable_regions,
        contrasts=None,
        adjpval_thr=0.05,
        log2fc_thr=np.log2(1.5),
        n_cpu=5,
        split_pattern='-'
    )

    # Remove empty markers
    for DAR in list(markers_dict.keys()):
        if len(markers_dict[DAR]) == 0:
            del markers_dict[DAR]

    # Save topics
    folder = os.path.join(par['temp_dir'], 'region_sets', 'Topics_otsu')
    os.makedirs(folder, exist_ok=True)
    for topic in binarized_cell_topic:
        region_names_to_coordinates(
            region_bin_topics_otsu[topic].index
        ).sort_values(
            ['Chromosome', 'Start', 'End']
        ).to_csv(
            os.path.join(folder, f'{topic}.bed'),
            sep='\t',
            header=False, index=False
        )
    folder = os.path.join(par['temp_dir'], 'region_sets', 'Topics_top_3k')
    os.makedirs(folder, exist_ok=True)
    for topic in region_bin_topics_top_3k:
        region_names_to_coordinates(
            region_bin_topics_top_3k[topic].index
        ).sort_values(
            ['Chromosome', 'Start', 'End']
        ).to_csv(
            os.path.join(folder, f'{topic}.bed'),
            sep='\t',
            header=False, index=False
        )
    folder = os.path.join(par['temp_dir'], 'candidate_enhancers')
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, 'region_bin_topics_otsu.pkl'), 'wb') as f:
        pickle.dump(region_bin_topics_otsu, f)
    with open(os.path.join(folder, 'region_bin_topics_top3k.pkl'), 'wb') as f:
        pickle.dump(region_bin_topics_top_3k, f)

    # Save DARs
    folder = os.path.join(par['temp_dir'], 'region_sets', 'DARs_cell_type')
    os.makedirs(folder, exist_ok=True)
    for cell_type in markers_dict:
        region_names_to_coordinates(
            markers_dict[cell_type].index
        ).sort_values(
            ['Chromosome', 'Start', 'End']
        ).to_csv(
            os.path.join(folder, f'{cell_type}.bed'),
            sep='\t',
            header=False,
            index=False
        )
    folder = os.path.join(par['temp_dir'], 'candidate_enhancers')
    with open(os.path.join(folder, 'markers_dict.pkl'), 'wb') as f:
        pickle.dump(markers_dict, f)

    # Get gene activity
    pr_annotation = pd.read_table(
        os.path.join(par['temp_dir'], 'qc', 'tss.bed')
    ).rename(
        {'Name': 'Gene', '# Chromosome': 'Chromosome'}, axis=1)
    pr_annotation['Transcription_Start_Site'] = pr_annotation['Start']
    pr_annotation = pr.PyRanges(pr_annotation)
    gene_act, weights = get_gene_activity(
        imputed_acc_obj,
        pr_annotation,
        chromsizes,
        use_gene_boundaries=True,
        upstream=[1000, 100000],
        downstream=[1000, 100000],
        distance_weight=True,
        decay_rate=1,
        extend_gene_body_upstream=10000,
        extend_gene_body_downstream=500,
        gene_size_weight=False,
        gene_size_scale_factor='median',
        remove_promoters=False,
        average_scores=True,
        scale_factor=1,
        extend_tss=[10, 10],
        gini_weight = True,
        return_weights= True
    )

    # Infer differentially accessible genes
    DAG_markers_dict = find_diff_features(
        cistopic_obj,
        gene_act,
        variable='cell_type',
        var_features=None,
        contrasts=None,
        adjpval_thr=0.05,
        log2fc_thr=np.log2(1.5),
        n_cpu=par['num_workers'],
        split_pattern='-'
    )

    # Exporting to loom
    os.makedirs(os.path.join(par['temp_dir'], 'loom'), exist_ok=True)
    cluster_markers = {'cell_type': markers_dict}
    filepath = os.path.join(par['temp_dir'], 'loom', f'region-accessibility.loom')
    print(f'Saving region accessibility to {filepath}')
    export_region_accessibility_to_loom(
        accessibility_matrix=imputed_acc_obj,
        cistopic_obj=cistopic_obj,
        binarized_topic_region=region_bin_topics_otsu,
        binarized_cell_topic=binarized_cell_topic,
        out_fname=filepath,
        cluster_annotation=['cell_type'],
        cluster_markers=cluster_markers,
        nomenclature='hg38',
        split_pattern='-'
    )
    filepath = os.path.join(par['temp_dir'], 'loom', f'gene-activity.loom')
    print(f'Saving gene activity to {filepath}')
    export_gene_activity_to_loom(
        gene_activity_matrix=gene_act,
        cistopic_obj=cistopic_obj,
        out_fname=filepath,
        cluster_annotation=['cell_type'],
        cluster_markers=cluster_markers,
        nomenclature='hg38',
        split_pattern='-'
    )

def download_databases(par):
    def download(url: str, filepath: str) -> None:
        if os.path.exists(filepath):
            return
        print(f'Download {url}...')
        urlretrieve(url, filepath)
    # Download list of blacklist regions
    print('Download list of blacklist regions')
    url = 'https://raw.githubusercontent.com/aertslab/pycisTopic/d6a2f8c832c14faae07def1d3ad8755531f50ad5/blacklist/hg38-blacklist.v2.bed'
    if not os.path.exists(par['blacklist_path']):
        # response = requests.get(url)
        # with open(par['blacklist_path'], 'w') as f:
        #     f.write(response.text)
        download(url, par['blacklist_path'])
    url = 'https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl'
    if not os.path.exists(par['motif_annotation']):
        download(url, par['motif_annotation'])
    if not os.path.exists(par['SCORES_DB_PATH']):
        # download scores
        url = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather'
        download(url, par['SCORES_DB_PATH'])
    if not os.path.exists(par['RANKINGS_DB_PATH']):
        # download ranking
        url = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
        download(url, par['RANKINGS_DB_PATH'])
    if False: # create custom databases
        if not (os.path.exists(RANKINGS_DB_PATH) and os.path.exists(SCORES_DB_PATH)):

            # Download create_cisTarget_databases
            os.makedirs(os.path.join(par['temp_dir'], 'cistarget-db'), exist_ok=True)
            if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'create_cisTarget_databases')):
                with contextlib.chdir(os.path.join(par['temp_dir'], 'cistarget-db')):
                    subprocess.run(['git', 'clone', 'https://github.com/aertslab/create_cisTarget_databases'])
            # Download cluster-buster
            if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'cbust')):
                urlretrieve('https://resources.aertslab.org/cistarget/programs/cbust', os.path.join(par['temp_dir'], 'cistarget-db', 'cbust'))
            subprocess.run(['chmod', 'a+x', os.path.join(par['temp_dir'], 'cistarget-db', 'cbust')])

            # Download motif collection
            if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'v10nr_clust_public')):
                urlretrieve(
                    'https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip',
                    os.path.join(par['temp_dir'], 'cistarget-db', 'v10nr_clust_public.zip')
                )
                with zipfile.ZipFile(os.path.join(par['temp_dir'], 'cistarget-db', 'v10nr_clust_public.zip'), 'r') as zip_ref:
                    zip_ref.extractall(os.path.join(par['temp_dir'], 'cistarget-db'))

            # Download chromosome sizes
            if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.chrom.sizes')):
                urlretrieve(
                    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
                    os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.chrom.sizes')
                )

            # Download reference genome
            if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.fa')):
                print('Downloading reference genome', flush=True)
                if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.fa.gz')):
                    urlretrieve(
                        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
                        os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.fa.gz')
                    )
                with gzip.open(os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.fa.gz'), 'rb') as f_in:
                    with open(os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.fa'), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

            # Prepare fasta from consensus regions
            if not os.path.exists(os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.with_1kb_bg_padding.fa')):
                subprocess.run([
                    os.path.join(par['temp_dir'], 'cistarget-db', 'create_cisTarget_databases', 'create_fasta_with_padded_bg_from_bed.sh'),
                    os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.fa'),
                    os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.chrom.sizes'),
                    os.path.join(par['temp_dir'], 'consensus_peak_calling', 'consensus_regions.bed'),
                    os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.with_1kb_bg_padding.fa'),
                    '1000',
                    'yes'
                ])
            # Create cistarget databases
            with open(os.path.join(par['temp_dir'], 'cistarget-db', 'motifs.txt'), 'w') as f:
                for filename in os.listdir(os.path.join(par['temp_dir'], 'cistarget-db', 'v10nr_clust_public', 'singletons')):
                    f.write(f'{filename}\n')
            print(os.path.join(par['temp_dir'], 'cistarget-db', 'create_cisTarget_databases', 'create_cistarget_motif_databases.py'))
            subprocess.run('file ', os.path.join(par['temp_dir'], 'cistarget-db', 'create_cisTarget_databases', 'create_cistarget_motif_databases.py'))
            # aa
            # with contextlib.chdir(os.path.join(par['temp_dir'], 'cistarget-db')):
            result = subprocess.run([
                'python',
                os.path.join(par['temp_dir'], 'cistarget-db', 'create_cisTarget_databases', 'create_cistarget_motif_databases.py'),
                '-f', os.path.join(par['temp_dir'], 'cistarget-db', 'hg38.with_1kb_bg_padding.fa'),
                '-M', os.path.join(par['temp_dir'], 'cistarget-db', 'v10nr_clust_public', 'singletons'),
                '-m', os.path.join(par['temp_dir'], 'cistarget-db', 'motifs.txt'),
                '-c', os.path.join(par['temp_dir'], 'cistarget-db', 'cbust'),
                '-o', 'db',
                '--bgpadding', '1000',
                '-t', str(par['num_workers'])
            ], capture_output=True, text=True)
            # Print the result for debugging
            print(result.stdout)
            print(result.stderr)
    
def preprocess_rna(par):
    os.makedirs(os.path.join(par['temp_dir'], 'scRNA'), exist_ok=True)
    print("Preprocess RNA-seq", flush=True)
    # Load scRNA-seq data
    print("Load scRNA-seq data")
    adata_rna = anndata.read_h5ad(par['rna'])
    # Only keep cells with at least 200 expressed genes, and genes with at least 3 cells expressing them
    sc.pp.filter_cells(adata_rna, min_genes=200)
    sc.pp.filter_genes(adata_rna, min_cells=3)

    # Filter out doublets using scrublet
    # print("Filter out doublets using scrublet")
    # sc.external.pp.scrublet(adata_rna)
    # adata_rna = adata_rna[adata_rna.obs['predicted_doublet'] == False]

    # Filter based on mitochondrial and total counts
    # adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
    # sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Normalize data but keep track of original data
    adata_rna.raw = adata_rna
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)
    # Change barcodes to match the barcodes in the scATAC-seq data
    bar_codes = [f'{obs_name.replace("-", "")}-{donor_id}' for obs_name, donor_id in zip(adata_rna.obs_names, adata_rna.obs.donor_id)]
    adata_rna.obs_names = bar_codes
    # Save scRNA-seq data
    adata_rna.write_h5ad(os.path.join(par['temp_dir'], 'rna.h5ad'))

def snakemake_pipeline(par):
    import os
    snakemake_dir = os.path.join(par['temp_dir'], 'scplus_pipeline', 'Snakemake')
    if os.path.exists(snakemake_dir):
        import shutil
        shutil.rmtree(snakemake_dir)  

    os.makedirs(os.path.join(par['temp_dir'], 'scplus_pipeline'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'scplus_pipeline', 'temp'), exist_ok=True)
    
    pipeline_dir = os.path.join(par['temp_dir'], 'scplus_pipeline')
    # if not os.path.exists(pipeline_dir):
    subprocess.run(['scenicplus', 'init_snakemake', '--out_dir', pipeline_dir])
    print('snake make initialized', flush=True)

    # Load pipeline settings
    with open(os.path.join(snakemake_dir, 'config', 'config.yaml'), 'r') as f:
        settings = yaml.safe_load(f)
    print('output_data:', settings['output_data'], flush=True)

    # Update settings
    cwd = os.getcwd()
    print(cwd)
    
    settings['input_data']['cisTopic_obj_fname'] = f"{cwd}/{par['cistopic_object']}"
    settings['input_data']['GEX_anndata_fname'] = f"{cwd}/{os.path.join(par['temp_dir'], 'rna.h5ad')}"
    settings['input_data']['region_set_folder'] = f"{cwd}/{os.path.join(par['temp_dir'], 'region_sets')}"
    settings['input_data']['ctx_db_fname'] = f"{cwd}/{par['RANKINGS_DB_PATH']}"
    settings['input_data']['dem_db_fname'] = f"{cwd}/{par['SCORES_DB_PATH']}"
    settings['input_data']['path_to_motif_annotations'] = f"{cwd}/{par['motif_annotation']}"
    settings['params_general']['temp_dir'] = f"{cwd}/{os.path.join(par['temp_dir'], 'scplus_pipeline', 'temp')}"
    settings['params_general']['n_cpu'] = par['num_workers']
    settings['params_inference']['quantile_thresholds_region_to_gene'] = '0.7 0.8 0.9'
    settings['params_inference']['top_n_regionTogenes_per_gene'] = '10 15 25'
    settings['params_inference']['region_to_gene_importance_method'] = 'GBM'
    settings['params_inference']['tf_to_gene_importance_method'] = 'GBM'
    settings['params_inference']['rho_threshold'] = 0.0
    settings['params_inference']['region_to_gene_correlation_method'] = 'SR'
    settings['params_inference']['min_target_genes'] = 10
    settings['params_motif_enrichment']['species'] = 'homo_sapiens'
    settings['params_motif_enrichment']['motif_similarity_fdr'] = 0.00001
    settings['params_motif_enrichment']['dem_adj_pval_thr'] = 0.05
    settings['params_motif_enrichment']['dem_log2fc_thr'] = 0.5
    settings['params_motif_enrichment']['dem_promoter_space'] = 1000
    settings['params_motif_enrichment']['dem_motif_hit_thr'] = 3.0
    settings['params_motif_enrichment']['dem_max_bg_regions'] = 3000
    settings['params_motif_enrichment']['ctx_auc_threshold'] = 0.005
    settings['params_motif_enrichment']['ctx_nes_threshold'] = 3.0
    settings['params_data_preparation']['bc_transform_func'] = '\"lambda x: f\'{x}\'\"'  # '"lambda x: f''{x}''"'
    settings['params_data_preparation']['is_multiome'] = True
    settings['params_data_preparation']['key_to_group_by'] = ''
    settings['params_data_preparation']['nr_cells_per_metacells'] = 5
    settings['params_data_preparation']['species'] = 'hsapiens'
    settings['output_data']['scplus_mdata'] = f"{cwd}/{par['scplus_mdata']}"
    
    import os

    # List of file paths to check
    file_paths = [
        settings['input_data']['cisTopic_obj_fname'],
        settings['input_data']['GEX_anndata_fname'],
        settings['input_data']['region_set_folder'],
        settings['input_data']['ctx_db_fname'],
        settings['input_data']['dem_db_fname'],
        settings['input_data']['path_to_motif_annotations']
    ]

    # Assert all files exist
    for file_path in file_paths:
        assert os.path.exists(file_path), f"File not found: {file_path}"
    

    # Save pipeline settings
    print('output_data:', settings['output_data'], flush=True)
    with open(os.path.join(snakemake_dir, 'config', 'config.yaml'), 'w') as f:
        yaml.dump(settings, f)

    # Run pipeline
    print('run snakemake ', flush=True)

    with contextlib.chdir(snakemake_dir):
        # this fails to download atumatically so we do manually
        if True:
            url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"
            response = requests.get(url)
            with open('hg38.chrom.sizes', 'wb') as file:
                file.write(response.content)
        
        subprocess.run([
            f'snakemake',
            'touch'
        ])
        
        subprocess.run([
            f'snakemake',
            '--unlock'
        ])

        subprocess.run([
            f'snakemake',
            '--cores', str(par['num_workers'])])

def post_process(par):
    scplus_mdata = mudata.read(par['scplus_mdata'])

    # Save results
    grn_extended = scplus_mdata.uns['direct_e_regulon_metadata']

    grn_extended = grn_extended[['TF', 'Gene', 'rho_TF2G', 'Region']].drop_duplicates().reset_index(drop=True)
    grn_extended = grn_extended[grn_extended.rho_TF2G!=0]

    grn_extended.columns = ['source', 'target', 'weight', 'peak']

    grn_extended = grn_extended[['source', 'target', 'weight', 'peak']].drop_duplicates(ignore_index=True)

    net = grn_extended.groupby(['source', 'target'], as_index=False)['weight'].sum()

    grn_extended.to_csv(par['grn_extended'])

    net = net.astype(str)
    return net