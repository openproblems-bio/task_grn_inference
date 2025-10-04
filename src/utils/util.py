import pandas as pd
import anndata as ad
import numpy as np
from tqdm import tqdm
import scipy.sparse as sp

def naming_convention(dataset, method): 
    if (dataset in ['replogle', 'parsescience', 'xaira_HEK293T']) & (method in ['scprint']):
        dataset = f'{dataset}_sc'
    return f'{dataset}.{method}.{method}.prediction.h5ad'

def binarize_weight(weight):
    if weight > 0:
        return 1
    elif weight < 0:
        return -1
    else:
        return 0

def manage_layer(adata, par):
    dataset = adata.uns['dataset_id']
    if 'layer' in par:
        layer = par['layer']
    else:
        if 'lognorm' in adata.layers:
            layer = 'lognorm'
        elif 'X_norm' in adata.layers:
            layer = 'X_norm'
        else:
            raise ValueError('No layer specified and no default layer found (lognorm or X_norm)')
    if layer not in adata.layers:
        if ('X_norm' in adata.layers) & (dataset in ['adamson', 'norman', 'nakatake']):
            layer = 'X_norm'
            print(f'Layer {par["layer"]} not found, using X_norm instead')
        else:
            raise ValueError(f'Layer {par["layer"]} not found in the data. Available layers: {list(adata.layers.keys())}')
    return layer
def read_prediction(par):
    adata = ad.read_h5ad(par['prediction'])
    net = adata.uns['prediction']
    processed_net = process_links(net, par)
    return processed_net

def format_save_score(output, method_id, dataset_id, score_file):
    
    print('Write output to file', flush=True)
    print(output)
    metric_ids = output.columns.to_numpy()
    metric_values = output.values[0]

    output = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": dataset_id,
            "method_id": method_id,
            "metric_ids": metric_ids,
            "metric_values": metric_values
        }
    )
    output.write_h5ad(score_file, compression='gzip')
    print('Completed', flush=True)
def parse_args(par):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--rna', type=str, help='Path to the input RNA data in h5ad format.')
    parser.add_argument('--rna_all', type=str, help='Path to the input RNA all data in h5ad format.')
    parser.add_argument('--atac', type=str, help='Path to the input ATAC data in h5ad format.')
    parser.add_argument('--prediction', type=str, help='Path to the output prediction in h5ad format.')
    parser.add_argument('--score', type=str, help='Path to the output score in h5ad format.')
    parser.add_argument('--layer', type=str)
    parser.add_argument('--temp_dir', type=str)
    parser.add_argument('--tf_all', type=str)
    parser.add_argument('--skeleton', type=str)
    parser.add_argument('--apply_skeleton', action='store_true')
    parser.add_argument('--apply_tf', action='store_true')
    parser.add_argument('--max_n_links', type=int)
    parser.add_argument('--reg_type', type=str)
    parser.add_argument('--num_workers', type=int)
    parser.add_argument('--regulators_consensus', type=str)
    parser.add_argument('--evaluation_data', type=str)
    parser.add_argument('--ws_consensus', type=str)
    parser.add_argument('--ws_distance_background', type=str)
    parser.add_argument('--group_specific', type=str)

   
    
    args = parser.parse_args()
    for k, v in vars(args).items():
        if v is not None:
            par[k] = v
    return par

def verbose_print(verbose_level, message, level):
    if level <= verbose_level:
        print(message)


def verbose_tqdm(iterable, desc, level, verbose_level):
    if level <= verbose_level:
        return tqdm(iterable, desc=desc)
    else:
        return iterable  # Return the iterable without a progress bar


def basic_qc(
    adata, min_genes_per_cell=200, max_genes_per_cell=5000, min_cells_per_gene=10
):
    mt = adata.var_names.str.startswith("MT-")
    print("shape before ", adata.shape)
    # 1. stats
    total_counts = adata.X.sum(axis=1)
    n_genes_by_counts = (adata.X > 0).sum(axis=1)
    # mt_frac = adata[:, mt].X.sum(axis=1) / total_counts

    low_gene_filter = n_genes_by_counts < min_genes_per_cell
    high_gene_filter = n_genes_by_counts > max_genes_per_cell
    # mt_filter = mt_frac > max_mt_frac

    # 2. Filter cells
    # print(f'Number of cells removed: below min gene {low_gene_filter.sum()}, exceed max gene {high_gene_filter.sum()}')
    mask_cells = (~low_gene_filter) & (~high_gene_filter)
    #  (~mt_filter)
    # 3. Filter genes
    n_cells = (adata.X != 0).sum(axis=0)
    mask_genes = n_cells > min_cells_per_gene
    adata_f = adata[mask_cells, mask_genes]
    print("shape after ", adata_f.shape)
    return adata_f


def process_links(net, par):
    import numpy as np
    print('Original net shape: ', net.shape)
    if 'cell_type' in net.columns:
        print('Prediction contains cell type specific links. ')
    net_org = net.copy()
    cols = ["source", "target", "weight"]
    supp_cols = []

    if 'group_specific' in par and par['group_specific'] is not None:
        if par['group_specific'] in net.columns:
            supp_cols.append(par['group_specific'])
            
    # Check for symmetric links
    def per_group_process(net):
        flipped = net.rename(columns={"source": "target", "target": "source"})
        merged = net.merge(flipped, on=cols, how="inner")
        if not merged.empty:
            print("Warning: The network contains at least one symmetric link.")
        # Remove duplicates
        # net = net.drop_duplicates(subset=cols)
        
        # Remove self-loops
        net["source"] = net["source"].astype(str)
        net["target"] = net["target"].astype(str)
        net = net[net["source"] != net["target"]]
        
        # Aggregate duplicates by mean weight
        net["weight"] = pd.to_numeric(net["weight"], errors='coerce')
        net = net.groupby(["source", "target"], as_index=False)["weight"].mean()
        
        print(f"Network shape after cleaning: {net.shape}")
        # Limit the number of links
        max_links = par.get("max_n_links", 50000)
        # sort by absolute weight descending
        net = net.reindex(net["weight"].abs().sort_values(ascending=False).index).head(max_links)
        print(f"Network shape applying max_n_links: {net.shape}")
        return net
    print('Supplementary columns for grouping:', supp_cols)
    if len(supp_cols) == 0:
        net = per_group_process(net)
    else:
        print('Processing group-specific networks for:', supp_cols)
        net_store = []
        for sup_col in supp_cols:
            print(f"Processing group: {sup_col}")
            net = net_org[net_org[par['group_specific']] == sup_col]
            net = per_group_process(net)
            net[par['group_specific']] = sup_col
            net_store.append(net)
        net = pd.concat(net_store)
        net = net.reset_index(drop=True)
        print(f"Final network shape after group-specific processing: {net.shape}")
    # print(net)
    return net


def efficient_melting(net, gene_names, tf_all=None, symmetric=True):
    """Efficiently converts a network matrix into a DataFrame. If symmetric, only the upper triangle is considered.
    If not symmetric, all nonzero values are considered and rows are treated as source.
    If tf_all is not None, only the interactions with source as TFs are kept.
    """
    if symmetric:
        upper_triangle_indices = np.triu_indices_from(net, k=1)
        sources = np.array(gene_names)[upper_triangle_indices[0]]
        targets = np.array(gene_names)[upper_triangle_indices[1]]
        weights = net[upper_triangle_indices]
    else:
        row_indices, col_indices = np.where(net != 0)  # Extract all nonzero values
        sources = np.array(gene_names)[row_indices]
        targets = np.array(gene_names)[col_indices]
        weights = net[row_indices, col_indices]
    if tf_all is not None:
        mask_tf = np.isin(sources, tf_all)
        sources = sources[mask_tf]
        targets = targets[mask_tf]
        weights = weights[mask_tf]

    data = np.column_stack((sources, targets, weights))
    net_df = pd.DataFrame(data, columns=["source", "target", "weight"])

    return net_df


def corr_net(adata, tf_all, par) -> pd.DataFrame:
    # - read data
    X = adata.layers[par['layer']]
    if hasattr(X, "todense"):
        X = X.todense().A
    # - remove genes with 0 standard deviation
    gene_std = np.std(X, axis=0)
    nonzero_std_genes = gene_std > 0
    X = X[:, nonzero_std_genes]
    gene_names = adata[:, nonzero_std_genes].var_names.to_numpy()
    # - calculate correlation
    corr_type = par.get("corr_type", "pearson")
    if corr_type == "pearson":
        net = np.corrcoef(X.T)
    elif corr_type == "spearman":
        from scipy.stats import spearmanr
        rho, _ = spearmanr(X, axis=0)  # (n_genes x n_genes)
        net = rho
    else:
        raise ValueError(f"Unknown corr_type: {corr_type}")
    # - melt the matrix
    net = efficient_melting(net, gene_names, symmetric=True)

    # - subset to known TFs
    if par["apply_tf_methods"]:
        tf_all = np.intersect1d(tf_all, gene_names)
        net = net[net["source"].isin(tf_all)]

    # - process links: size control
    net = process_links(net, par)
    net = net.reset_index(drop=True)

    return net


# def plot_heatmap(scores, ax=None, name="", fmt="0.02f", cmap="viridis"):
#     import matplotlib.pyplot as plt
#     import numpy as np
#     import seaborn

#     if ax is None:
#         fig, ax = plt.subplots(1, 1, figsize=(4, 4), sharey=True)

#     # Normalize each column individually
#     scores_normalized = scores.apply(
#         lambda x: (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)), axis=0
#     )
#     scores_normalized = scores_normalized.round(2)
#     # scores_normalized['Rank'] = scores['Rank'].max()-scores['Rank']
#     # scores_normalized['Rank'] = scores_normalized['Rank']/scores_normalized['Rank'].max()

#     # Define the color scale range for each column (0 to 1 after normalization)
#     vmin = 0
#     vmax = 1

#     # Plot the heatmap with normalized scores
#     seaborn.heatmap(
#         scores_normalized,
#         ax=ax,
#         square=False,
#         cbar=False,
#         annot=True,
#         fmt=fmt,
#         vmin=vmin,
#         vmax=vmax,
#         cmap=cmap,
#     )
#     # Overlay the original (unnormalized) scores as annotations
#     # scores['Rank'] = scores['Rank'].astype(int)
#     # print(scores['Rank'])
#     # Overlay the original (unnormalized) scores as annotations
#     for text, (i, j) in zip(ax.texts, np.ndindex(scores.shape)):
#         value = scores.iloc[i, j]
#         if isinstance(value, np.int64):  # Check if the value is an integer for 'Rank'
#             text.set_text(f"{value:d}")
#         else:
#             text.set_text(f"{value:.2f}")

#     # Customize the axes and title
#     ax.tick_params(left=False, bottom=False)
#     ax.xaxis.set_tick_params(width=0)
#     ax.yaxis.set_tick_params(width=0)
#     ax.set_title(name, pad=10)

#     ax.xaxis.set_label_position("top")
#     ax.xaxis.tick_top()
#     ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="left")


def read_gmt(file_path: str) -> dict[str, list[str]]:
    """Reas gmt file and returns a dict of gene"""
    gene_sets = {}
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            gene_set_name = parts[0]
            gene_set_description = parts[1]
            genes = parts[2:]
            gene_sets[gene_set_name] = {
                "description": gene_set_description,
                "genes": genes,
            }
    return gene_sets


def read_gene_annotation(annotation_file):
    """
    Read the gene annotation file and extract TSS
    """

    gtf_columns = [
        "chromosome",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]
    gtf_df = pd.read_csv(
        annotation_file, sep="\t", comment="#", names=gtf_columns, low_memory=False
    )

    gtf_df["gene_name"] = gtf_df["attributes"].str.extract(r'gene_name "([^"]+)"')

    # Keep only gene-level annotations
    annotation_df = gtf_df[gtf_df["feature"] == "gene"][
        ["chromosome", "start", "end", "strand", "gene_name"]
    ]

    # Compute TSS
    annotation_df["TSS"] = annotation_df.apply(
        lambda x: x["start"] if x["strand"] == "+" else x["end"], axis=1
    )

    return annotation_df[["chromosome", "start", "end", "TSS", "strand", "gene_name"]]




def fetch_gene_info():
    from pybiomart import Server

    # Connect to Ensembl server
    server = Server(host="http://www.ensembl.org")

    # Select the dataset for human genes
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

    # Query relevant attributes
    df = dataset.query(
        attributes=[
            "ensembl_gene_id",  # Ensembl Gene ID
            "external_gene_name",  # Gene Name
            "chromosome_name",  # Chromosome
            "start_position",  # Start site
            "end_position",  # End site
        ]
    )
    df.columns = ["gene_id", "gene_name", "chr", "start", "end"]
    # df = df[df['chr'].isin(['X', 'Y']+list(map(str, range(1, 23))))]
    # - keep those genes with one mapping
    df = df.groupby("gene_name").filter(lambda x: len(x) == 1)
    df.set_index("gene_name", inplace=True)

    return df

def download_annotation(par):
    import os
    import requests

    if not os.path.exists(par['annotation_file']):
        print("Downloading prior started")
        response = requests.get("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz")
        if response.status_code == 200:
            with open(par['annotation_file'], 'wb') as file:
                file.write(response.content)
            print(f"File downloaded and saved to {par['annotation_file']}")
        else:
            print(f"Failed to download the gencode.v45.annotation.gtf.gz. Status code: {response.status_code}")
        print("Downloading prior ended")

def read_gtf_as_df(gtf_path: str) -> pd.DataFrame:
    """
    Read a GTF/GFF3 file (plain or gzipped) and return gene-level annotation
    with TSS included.

    Returns
    -------
    pd.DataFrame
        Columns: ['Chromosome', 'Start', 'End', 'Strand', 'Gene', 'Transcription_Start_Site']
    """
    import gzip
    import pandas as pd

    rows = []

    # Open file with automatic gzip detection
    if gtf_path.endswith(".gz"):
        f = gzip.open(gtf_path, "rt")
    else:
        try:
            f = open(gtf_path, "r", encoding="utf-8")
            f.readline()
            f.seek(0)
        except UnicodeDecodeError:
            f = gzip.open(gtf_path, "rt")

    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if parts[2] != "gene":
            continue
        chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]

        # attributes
        attr_str = parts[8]
        attrs = {}
        for attr in attr_str.strip().split(";"):
            if attr.strip() == "":
                continue
            key, value = attr.strip().split(" ", 1)
            attrs[key] = value.strip('"')

        gene_name = attrs.get("gene_name", attrs.get("gene_id", ""))

        # Compute TSS
        tss = start if strand == "+" else end

        rows.append([chrom, start, end, strand, gene_name, tss])

    f.close()

    df = pd.DataFrame(
        rows,
        columns=["Chromosome", "Start", "End", "Strand", "Gene", "Transcription_Start_Site"]
    )
    return df