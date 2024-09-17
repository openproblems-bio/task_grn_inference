import os
import sys
import yaml
import pickle
import tempfile
import contextlib
import hashlib
import shutil
import requests
import traceback
import subprocess
import gc
import gzip
import tarfile
from urllib.request import urlretrieve

import numpy as np
import scanpy as sc
import pandas as pd
import anndata
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget
import polars
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


## VIASH START
par = {
  'multiomics_rna': 'resources_test/grn-benchmark/multiomics_rna.h5ad',
  'multiomics_atac': 'resources_test/grn-benchmark/multiomics_atac.h5ad',
  'temp_dir': 'output/scenicplus',
  'cistopic_object': 'output/scenicplus/pycistopic.pkl',
  'prediction': 'output/scenicplus/prediction.csv',
  'qc': False,
  'num_workers': 20,
  'cell_topic': 'output/scenicplus/cell_topic.csv',
}
## VIASH END
## installing missing software. #TODO: fix this later
try:
    subprocess.run(['yum', '-y', 'install', 'htslib'], check=True)
    print("htslib installed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Installation failed: {e}")

# Bug in pycistopic: Import is missing in pycisTopic.loom,
# so TSNE must be dynamically added to the library's namespace.
setattr(pycisTopic.loom, 'TSNE', TSNE)

os.environ['MALLET_MEMORY'] = '200G'
try:
    sys.path.append(meta['resources_dir'])
except NameError:
    pass

par['temp_dir'] = par['temp_dir']
atac_dir = os.path.join(par['temp_dir'], 'atac')
os.makedirs(atac_dir, exist_ok=True)

# Get list of samples (e.g., donors)
print('Collect list of samples')
adata_atac = anndata.read_h5ad(par['multiomics_atac'])
unique_donor_ids = [s.replace(' ', '_') for s in adata_atac.obs.donor_id.cat.categories]
unique_cell_types = [s.replace(' ', '_') for s in adata_atac.obs.cell_type.cat.categories]
del adata_atac

def process_atac(par):
    print("---Run pre-process started ---", flush=True)
    # Create one individual ATAC-seq file per donor
    fragments_dict = {}
    for donor_id in unique_donor_ids:
        filepath = os.path.join(atac_dir, f'{donor_id}.tsv')
        if not os.path.exists(filepath):
            print(f'Create tsv file {filepath}')
            adata_atac = anndata.read_h5ad(par['multiomics_atac'])
            adata_atac.obs.cell_type = [s.replace(' ', '_') for s in adata_atac.obs.cell_type]
            adata_atac.obs.donor_id = [s.replace(' ', '_') for s in adata_atac.obs.donor_id]
            adata_atac_one_donor = adata_atac[adata_atac.obs.donor_id == donor_id]
            coo_matrix = adata_atac_one_donor.X.tocoo()
            rows = coo_matrix.row
            cols = coo_matrix.col
            counts = coo_matrix.data
            row_names = adata_atac_one_donor.obs.index[rows]
            coord_names = adata_atac_one_donor.var.index[cols]
            del adata_atac, adata_atac_one_donor
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
        if not os.path.exists(compressed_filepath):

            print(f'Sort and compress tsv file {filepath}')

            # Sort file by genomic coordinates
            sorted_filepath = filepath + '.sorted.tsv'
            os.system(f'sort -k1,1 -k2,2n {filepath} > {sorted_filepath}')

            # Compression
            subprocess.run(['bgzip', sorted_filepath, '-o', compressed_filepath])

        fragments_dict[donor_id] = compressed_filepath

        # Index file using tabix
        if not os.path.exists(compressed_filepath + '.tbi'):
            print(f'Index compressed file {compressed_filepath} using tabix')
            subprocess.run(['tabix', compressed_filepath, '-p', 'bed'])

    # Collect cell metadata
    print(f'Collect cell metadata')
    adata_atac = anndata.read_h5ad(par['multiomics_atac'])
    donor_ids = [s.replace(' ', '_') for s in adata_atac.obs.donor_id]
    index = [barcode.replace('-', '') + '-' + donor_id for donor_id, barcode in zip(donor_ids, adata_atac.obs.index)]
    cell_data = pd.DataFrame({
        'cell_type': [s.replace(' ', '_') for s in adata_atac.obs.cell_type.to_numpy()],
        'donor_id': [s.replace(' ', '_') for s in adata_atac.obs.donor_id.to_numpy()]
    }, index=index, copy=False)

    # Get chromosome sizes
    print(f'Download chromosome sizes from UCSC')
    chromsizes = pd.read_table(
        'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
        header=None,
        names=['Chromosome', 'End']
    )
    chromsizes.insert(1, 'Start', 0)

    print(f'Start pseudobulking')
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bed_files'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bw_files'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/MACS'), exist_ok=True)
    os.makedirs(os.path.join(par['temp_dir'], 'consensus_peak_calling/tmp'), exist_ok=True)

    # Check if individual fragment files are available for each cell type
    valid = True
    bed_paths = {}
    for cell_type in unique_cell_types:
        filepath = os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bed_files', f'{cell_type}.fragments.tsv.gz')
        bed_paths[cell_type] = filepath
        if not os.path.exists(filepath):
            valid = False

    # Perform pseudobulking for each cell type
    if not valid:
        print('Pseudobulk each cell type')
        _, bed_paths = export_pseudobulk(
            input_data=cell_data,
            variable='cell_type',
            sample_id_col='donor_id',
            chromsizes=chromsizes,
            bed_path=os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bed_files'),
            bigwig_path=os.path.join(par['temp_dir'], 'consensus_peak_calling/pseudobulk_bw_files'),
            path_to_fragments=fragments_dict,
            n_cpu=10,
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

    # Download list of blacklist regions
    if not os.path.exists(os.path.join(par['temp_dir'], 'hg38-blacklist.v2.bed')):
        print('Download list of blacklist regions')
        url = 'https://raw.githubusercontent.com/aertslab/pycisTopic/d6a2f8c832c14faae07def1d3ad8755531f50ad5/blacklist/hg38-blacklist.v2.bed'
        response = requests.get(url)
        with open(os.path.join(par['temp_dir'], 'hg38-blacklist.v2.bed'), 'w') as f:
            f.write(response.text)

    # Consensus peak calling
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peak_dict,
        peak_half_width=250,
        chromsizes=chromsizes,
        path_to_blacklist=os.path.join(par['temp_dir'], 'hg38-blacklist.v2.bed')
    )
    consensus_peaks.to_bed(
        path=os.path.join(par['temp_dir'], 'consensus_peak_calling/consensus_regions.bed'),
        keep=True,
        compression='infer',
        chain=False
    )

    # Download TSS annotations
    os.makedirs(os.path.join(par['temp_dir'], 'qc'), exist_ok=True)
    if not os.path.exists(os.path.join(par['temp_dir'], 'qc', 'tss.bed')):
        subprocess.run([
            'pycistopic', 'tss', 'get_tss',
            '--output', os.path.join(par['temp_dir'], 'qc', 'tss.bed'),
            '--name', 'hsapiens_gene_ensembl',
            '--to-chrom-source', 'ucsc',
            '--ucsc', 'hg38'
        ])
    print("---Run pre-process completed ---", flush=True)

def run_cistopic(par):
    if os.path.exists(os.path.join(par['temp_dir'], 'cistopic_obj.pkl')):
        print('Cistopic object exists. Skipping run_cistopic.')
        return 
    print("---Run cistopic---", flush=True)
    if par['qc']:  # Whether to perform quality control

        # Compute QC metrics
        print('Perform QC')
        for donor_id in unique_donor_ids:
            filepath = os.path.join(atac_dir, f'{donor_id}.tsv.gz')
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
                '--threads', '10'
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
                path_to_blacklist=os.path.join(par['temp_dir'], 'hg38-blacklist.v2.bed'),
                metrics=sample_metrics,
                valid_bc=sample_id_to_barcodes_passing_filters[sample_id],
                n_cpu=10,
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
                path_to_blacklist=os.path.join(par['temp_dir'], 'hg38-blacklist.v2.bed'),
                n_cpu=10,
                project=donor_id,
                split_pattern='-'
            )
            cistopic_obj_list.append(cistopic_obj)

    # Add metadata to cistopic objects
    for i in range(len(cistopic_obj_list)):
        cistopic_obj_list[i].add_cell_data(cell_data, split_pattern='-')

    # Infer doublets using scrublet
    print('Infer doublets using scrublet')
    for i in range(len(cistopic_obj_list)):
        scrub = scr.Scrublet(cistopic_obj_list[i].fragment_matrix.T, expected_doublet_rate=0.1)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        scrub.plot_histogram();
        scrub.call_doublets(threshold=0.22)
        scrub.plot_histogram();
        scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T
        cistopic_obj_list[i].add_cell_data(scrublet, split_pattern='-')

    # Combine samples into a single cistopic object
    if len(cistopic_obj_list) == 1:
        cistopic_obj = cistopic_obj_list[0]
    else:
        cistopic_obj = merge(cistopic_obj_list, is_acc=1, copy=False, split_pattern='-')
    
    # Save cistopic objects
    with open(os.path.join(par['temp_dir'], 'cistopic_obj.pkl'), 'wb') as f:
        pickle.dump(cistopic_obj, f)
    # save cell_topic csv 
    cell_topic = cistopic_obj.selected_model.cell_topic.T
    cell_topic.index = cell_topic.index.str.split('-').str[0]
    cell_topic.to_csv(par['cell_topic'])
    
    print("---Run cistopic completed---", flush=True)
     
def run_post_cistopic(par):
    print("---Run post cistopic---", flush=True)
    # Load cistopic objects
    with open(os.path.join(par['temp_dir'], 'cistopic_obj.pkl'), 'rb') as f:
        cistopic_obj = pickle.load(f)
        
    chromsizes = pd.read_table(os.path.join(par['temp_dir'], "qc", "hg38.chrom_sizes_and_alias.tsv"))
    chromsizes.rename({"# ucsc": "Chromosome", "length": "End"}, axis = 1, inplace = True)
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
    print("---Run post cistopic completed---", flush=True)


process_atac(par)
run_cistopic(par)
run_post_cistopic(par)