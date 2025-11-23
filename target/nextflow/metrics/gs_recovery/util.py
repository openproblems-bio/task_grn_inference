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
    pd.set_option('display.max_columns', None)
    print('Write output to file', flush=True)

    if not isinstance(output, pd.DataFrame):
        raise ValueError("Expected 'output' to be a pandas DataFrame.")

    metric_ids = output.columns.tolist()
    metric_values = output.to_numpy().astype(str)[0]  # assuming single row of metrics

    

    adata = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": dataset_id,
            "method_id": method_id,
            "metric_ids": metric_ids,
            "metric_values": metric_values  # store as numpy array directly
        }
    )

    adata.write_h5ad(score_file, compression='gzip')
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
    parser.add_argument('--evaluation_data_de', type=str)
    parser.add_argument('--evaluation_data_sc', type=str)
    parser.add_argument('--ground_truth_unibind', type=str)
    parser.add_argument('--ground_truth_chipatlas', type=str)
    parser.add_argument('--ground_truth_remap', type=str)
    parser.add_argument('--geneset_hallmark_2020', type=str, help='Hallmark 2020 geneset GMT file')
    parser.add_argument('--geneset_kegg_2021', type=str, help='KEGG 2021 geneset GMT file')
    parser.add_argument('--geneset_reactome_2022', type=str, help='Reactome 2022 geneset GMT file')
    parser.add_argument('--geneset_go_bp_2023', type=str, help='GO BP 2023 geneset GMT file')
    parser.add_argument('--geneset_bioplanet_2019', type=str, help='BioPlanet 2019 geneset GMT file')
    parser.add_argument('--geneset_wikipathways_2019', type=str, help='WikiPathways 2019 geneset GMT file')

   
    
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



def add_gene_id(adata):
    df_annot = fetch_gene_info()[['gene_id']]
    genes_1 = df_annot.index 
    genes_2 = adata.var.index
    common_genes = genes_1.intersection(genes_2)
    assert len(common_genes) > 0, "No common genes between annotation and adata.var"
    merged = adata.var.merge(df_annot, left_index=True, right_index=True, how='left')
    valid_genes = merged.dropna().index
    adata = adata[:, valid_genes].copy()
    adata.var = merged.loc[valid_genes]
    return adata
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


from typing import List

import numpy as np


class WeightedDegree(object):

    def __init__(self, node_idx: int, degree: int, weight: float):
        """Class used to sort nodes by in-degree or out-degree.

        Ties are broken by using the weights from the GRN.

        Args:
            node_idx
            degree: in-degree or out-degree of a node.
            weight: sum of weights of incoming or outcoming edges.
        """
        self.node_idx: int = int(node_idx)
        self.degree: int = int(degree)
        self.weight: float = float(weight)

    def __lt__(self, other: "WeightedDegree") -> bool:
        if self.degree < other.degree:
            return True
        elif self.degree == other.degree:
            return self.weight < other.weight
        else:
            return False

    def __eq__(self, other: "WeightedDegree") -> bool:
        return (self.degree == other.degree) and (self.weight == other.weight)

    @staticmethod
    def from_grn(A: np.ndarray, incoming: bool = True) -> List["WeightedDegree"]:
        if not incoming:
            A = A.T
        D = (A != 0)
        degrees = np.sum(D, axis=0)
        weights = np.sum(A, axis=0)
        return [WeightedDegree(i, degrees[i], weights[i]) for i in range(len(degrees))]

def create_grn_baseline(A):
    """
    Optimized deterministic baseline for a directed simple graph.
    Preserves in/out degree sequences of A (A may be weighted/signed; nonzeros indicate edges).
    Returns a weighted B with degree-preserving shuffled topology and reassigned weights.
    
    Optimizations:
    - Vectorized degree calculations
    - Pre-allocated arrays
    - Reduced nested loops
    - Early termination conditions
    """
    # Ensure no self-regulation
    A = np.copy(A)
    np.fill_diagonal(A, 0)
    n = A.shape[0]
    
    # Vectorized degree calculations
    A_binary = (A != 0).astype(np.int32)
    in_deg = A_binary.sum(axis=0)
    out_deg = A_binary.sum(axis=1)
    in_wgt = np.abs(A).sum(axis=0)
    out_wgt = np.abs(A).sum(axis=1)
    
    # Create sorting keys for deterministic ordering
    # Primary: degree (desc), Secondary: weight (desc), Tertiary: index (asc)
    in_order = np.lexsort((np.arange(n), -in_wgt, -in_deg))
    out_order = np.lexsort((np.arange(n), -out_wgt, -out_deg))
    
    # Pre-allocate baseline structure
    B = np.zeros_like(A)
    residual_in_deg = in_deg.copy()
    
    # Build topology (edges without weights)
    for u in out_order:
        k = out_deg[u]
        if k == 0:
            continue
            
        # Vectorized selection of valid targets
        valid_mask = (residual_in_deg > 0) & (np.arange(n) != u) & (B[u, :] == 0)
        valid_targets = np.where(valid_mask)[0]
        
        if len(valid_targets) < k:
            # Simple fallback: take what we can get
            k = len(valid_targets)
        
        if k > 0:
            # Sort by residual in-degree (prefer high-degree nodes)
            target_order = np.argsort(-residual_in_deg[valid_targets])[:k]
            picks = valid_targets[target_order]
            
            # Place edges
            B[u, picks] = 1
            residual_in_deg[picks] -= 1
    
    np.fill_diagonal(B, 0)
    
    # Assign weights deterministically
    # Cache initial degrees for target ranking
    init_in_deg = in_deg.copy()
    init_in_wgt = in_wgt.copy()
    
    # Process each source node
    for u in range(n):
        # Get original weights from A
        orig_mask = (A[u, :] != 0)
        if not orig_mask.any():
            continue
            
        orig_targets = np.where(orig_mask)[0]
        W = A[u, orig_targets]
        
        # Sort weights by |magnitude| desc, then value, then index
        order_w = np.lexsort((orig_targets, W, -np.abs(W)))
        W_sorted = W[order_w]
        
        # Get new targets from B
        new_targets = np.where(B[u, :] == 1)[0]
        if len(new_targets) == 0:
            continue
            
        # Rank new targets by original graph importance
        target_keys = np.lexsort((new_targets, -init_in_wgt[new_targets], -init_in_deg[new_targets]))
        targets_ranked = new_targets[target_keys]
        
        # Match counts
        n_weights = min(len(W_sorted), len(targets_ranked))
        if n_weights > 0:
            B[u, targets_ranked[:n_weights]] = W_sorted[:n_weights]
    
    return B

