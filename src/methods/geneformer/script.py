import os
import shutil
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from geneformer import TranscriptomeTokenizer
from huggingface_hub import hf_hub_download
from scipy.sparse import csr_matrix

## VIASH START
[par] = {
    "rna": "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
    "tf_all": "resources_test/prior/tf_all.csv",
    "prediction": "output/geneformer/prediction.h5ad",
    "max_n_links": 50000,
    "batch_size": 16,
    "temp_dir": "output/geneformer",
    "file_dir": "src/methods/geneformer",
    "num_genes": 5000,
    "max_cells": 2000,
}
## VIASH END
meta = {"name": "geneformer"}

# imports
import logging
import pickle
from collections import Counter
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from geneformer import perturber_utils as pu
from tdigest import TDigest
from tqdm.auto import trange

logger = logging.getLogger(__name__)


# extract embeddings
def get_embs(
    model,
    filtered_input_data,
    emb_mode,
    layer_to_quant,
    pad_token_id,
    forward_batch_size,
    token_gene_dict,
    special_token=False,
    summary_stat=None,
    silent=False,
    get_avg_attentions=False,
    genelist=None,
):
    model_input_size = pu.get_model_input_size(model)
    total_batch_length = len(filtered_input_data)

    if summary_stat is None:
        embs_list = []
    elif summary_stat is not None:
        # get # of emb dims
        emb_dims = pu.get_model_emb_dims(model)
        if emb_mode == "cell":
            # initiate tdigests for # of emb dims
            embs_tdigests = [TDigest() for _ in range(emb_dims)]
        if emb_mode == "gene":
            gene_set = list(
                {
                    element
                    for sublist in filtered_input_data["input_ids"]
                    for element in sublist
                }
            )
            # initiate dict with genes as keys and tdigests for # of emb dims as values
            embs_tdigests_dict = {
                k: [TDigest() for _ in range(emb_dims)] for k in gene_set
            }

    # Check if CLS and EOS token is present in the token dictionary
    cls_present = any("<cls>" in value for value in token_gene_dict.values())
    eos_present = any("<eos>" in value for value in token_gene_dict.values())
    gene_token_dict = {v: k for k, v in token_gene_dict.items()}
    if emb_mode == "cls":
        assert cls_present, "<cls> token missing in token dictionary"
        # Check to make sure that the first token of the filtered input data is cls token
        cls_token_id = gene_token_dict["<cls>"]
        assert filtered_input_data["input_ids"][0][0] == cls_token_id, (
            "First token is not <cls> token value"
        )
    elif emb_mode == "cell":
        if cls_present:
            logger.warning(
                "CLS token present in token dictionary, excluding from average."
            )
        if eos_present:
            logger.warning(
                "EOS token present in token dictionary, excluding from average."
            )

    overall_max_len = 0
    if get_avg_attentions:
        avg_attentions = torch.zeros((len(genelist), len(genelist)), device="cuda")
        attention_counts = torch.zeros_like(avg_attentions)

    for i in trange(0, total_batch_length, forward_batch_size, leave=(not silent)):
        max_range = min(i + forward_batch_size, total_batch_length)

        minibatch = filtered_input_data.select([i for i in range(i, max_range)])

        max_len = int(max(minibatch["length"]))
        original_lens = torch.tensor(minibatch["length"], device="cuda")
        minibatch.set_format(type="torch")

        input_data_minibatch = minibatch["input_ids"]
        input_data_minibatch = pu.pad_tensor_list(
            input_data_minibatch, max_len, pad_token_id, model_input_size
        )
        with torch.no_grad():
            outputs = model(
                input_ids=input_data_minibatch.to("cuda"),
                attention_mask=pu.gen_attention_mask(minibatch),
                output_attentions=get_avg_attentions,
            )

        embs_i = outputs.hidden_states[layer_to_quant]
        attentions = outputs.attentions
        # Calculate the average of all attention heads
        if get_avg_attentions:
            # avg across heads
            # sum across cells
            attentions = torch.mean(torch.stack(attentions), dim=(0, 2))
            attentions = attentions[:, 1:-1, 1:-1]
            del outputs
            torch.cuda.empty_cache()
            # Create a mapping from token IDs to genelist indices

            token_to_loc = {
                gene_token_dict[tok_id]: i
                for i, tok_id in enumerate(genelist)
                if tok_id in gene_token_dict
            }

            # Create a tensor of genelist indices for each token in the batch
            gene_indices = torch.tensor(
                [
                    [token_to_loc.get(token.item(), -1) for token in seq[1:-1]]
                    for seq in input_data_minibatch
                ],
                device=attentions.device,
            )

            # Create a mask for valid gene indices
            valid_mask = gene_indices != -1

            # Use advanced indexing to map attention values to gene pairs
            batch, seq_len, _ = attentions.shape
            row_indices = gene_indices.unsqueeze(-1).expand(-1, -1, seq_len)
            col_indices = gene_indices.unsqueeze(1).expand(-1, seq_len, -1)
            avg_attentions.index_put_(
                (row_indices[valid_mask], col_indices[valid_mask]),
                attentions[valid_mask],
                accumulate=True,
            )
        # Compute which values were filled to calculate the mean
        attention_counts.index_put_(
            (row_indices[valid_mask], col_indices[valid_mask]),
            torch.ones_like(attentions[valid_mask]),
            accumulate=True,
        )

        # This represents the average attention across all heads for each token pair

        if emb_mode == "cell":
            if cls_present:
                non_cls_embs = embs_i[:, 1:, :]  # Get all layers except the embs
                if eos_present:
                    mean_embs = pu.mean_nonpadding_embs(non_cls_embs, original_lens - 2)
                else:
                    mean_embs = pu.mean_nonpadding_embs(non_cls_embs, original_lens - 1)
            else:
                mean_embs = pu.mean_nonpadding_embs(embs_i, original_lens)
            if summary_stat is None:
                embs_list.append(mean_embs)
            elif summary_stat is not None:
                # update tdigests with current batch for each emb dim
                accumulate_tdigests(embs_tdigests, mean_embs, emb_dims)
            del mean_embs
        elif emb_mode == "gene":
            if summary_stat is None:
                embs_list.append(embs_i)
            elif summary_stat is not None:
                for h in trange(len(minibatch)):
                    length_h = minibatch[h]["length"]
                    input_ids_h = minibatch[h]["input_ids"][0:length_h]

                    # double check dimensions before unsqueezing
                    embs_i_dim = embs_i.dim()
                    if embs_i_dim != 3:
                        logger.error(
                            f"Embedding tensor should have 3 dimensions, not {embs_i_dim}"
                        )
                        raise

                    embs_h = embs_i[h, :, :].unsqueeze(dim=1)
                    dict_h = dict(zip(input_ids_h, embs_h))
                    for k in dict_h.keys():
                        accumulate_tdigests(
                            embs_tdigests_dict[int(k)], dict_h[k], emb_dims
                        )
                    del embs_h
                    del dict_h
        elif emb_mode == "cls":
            cls_embs = embs_i[:, 0, :].clone().detach()  # CLS token layer
            embs_list.append(cls_embs)
            del cls_embs

        overall_max_len = max(overall_max_len, max_len)
        if not get_avg_attentions:
            del outputs
        del minibatch
        del input_data_minibatch
        del embs_i

        torch.cuda.empty_cache()

    if summary_stat is None:
        if (emb_mode == "cell") or (emb_mode == "cls"):
            embs_stack = torch.cat(embs_list, dim=0)
        elif emb_mode == "gene":
            embs_stack = pu.pad_tensor_list(
                embs_list,
                overall_max_len,
                pad_token_id,
                model_input_size,
                1,
                pu.pad_3d_tensor,
            )

    # calculate summary stat embs from approximated tdigests
    elif summary_stat is not None:
        if emb_mode == "cell":
            if summary_stat == "mean":
                summary_emb_list = tdigest_mean(embs_tdigests, emb_dims)
            elif summary_stat == "median":
                summary_emb_list = tdigest_median(embs_tdigests, emb_dims)
            embs_stack = torch.tensor(summary_emb_list)
        elif emb_mode == "gene":
            if summary_stat == "mean":
                [
                    update_tdigest_dict_mean(embs_tdigests_dict, gene, emb_dims)
                    for gene in embs_tdigests_dict.keys()
                ]
            elif summary_stat == "median":
                [
                    update_tdigest_dict_median(embs_tdigests_dict, gene, emb_dims)
                    for gene in embs_tdigests_dict.keys()
                ]
            return embs_tdigests_dict
    if get_avg_attentions:
        # Avoid division by zero
        attention_counts = torch.clamp(attention_counts, min=1)
        # Calculate the mean attention
        avg_attentions = avg_attentions / attention_counts
    else:
        avg_attentions = None
    return embs_stack, avg_attentions


def accumulate_tdigests(embs_tdigests, mean_embs, emb_dims):
    # note: tdigest batch update known to be slow so updating serially
    [
        embs_tdigests[j].update(mean_embs[i, j].item())
        for i in range(mean_embs.size(0))
        for j in range(emb_dims)
    ]


def update_tdigest_dict(embs_tdigests_dict, gene, gene_embs, emb_dims):
    embs_tdigests_dict[gene] = accumulate_tdigests(
        embs_tdigests_dict[gene], gene_embs, emb_dims
    )


def update_tdigest_dict_mean(embs_tdigests_dict, gene, emb_dims):
    embs_tdigests_dict[gene] = tdigest_mean(embs_tdigests_dict[gene], emb_dims)


def update_tdigest_dict_median(embs_tdigests_dict, gene, emb_dims):
    embs_tdigests_dict[gene] = tdigest_median(embs_tdigests_dict[gene], emb_dims)


def summarize_gene_embs(h, minibatch, embs_i, embs_tdigests_dict, emb_dims):
    length_h = minibatch[h]["length"]
    input_ids_h = minibatch[h]["input_ids"][0:length_h]
    embs_h = embs_i[h, :, :].unsqueeze(dim=1)
    dict_h = dict(zip(input_ids_h, embs_h))
    [
        update_tdigest_dict(embs_tdigests_dict, k, dict_h[k], emb_dims)
        for k in dict_h.keys()
    ]


def tdigest_mean(embs_tdigests, emb_dims):
    return [embs_tdigests[i].trimmed_mean(0, 100) for i in range(emb_dims)]


def tdigest_median(embs_tdigests, emb_dims):
    return [embs_tdigests[i].percentile(50) for i in range(emb_dims)]


def label_cell_embs(embs, downsampled_data, emb_labels):
    embs_df = pd.DataFrame(embs.cpu().numpy())
    if emb_labels is not None:
        for label in emb_labels:
            emb_label = downsampled_data[label]
            embs_df[label] = emb_label
    return embs_df


def label_gene_embs(embs, downsampled_data, token_gene_dict):
    gene_set = {
        element for sublist in downsampled_data["input_ids"] for element in sublist
    }
    gene_emb_dict = {k: [] for k in gene_set}
    for i in range(embs.size()[0]):
        length = downsampled_data[i]["length"]
        dict_i = dict(
            zip(
                downsampled_data[i]["input_ids"][0:length],
                embs[i, :, :].unsqueeze(dim=1),
            )
        )
        for k in dict_i.keys():
            gene_emb_dict[k].append(dict_i[k])
    for k in gene_emb_dict.keys():
        gene_emb_dict[k] = (
            torch.squeeze(torch.mean(torch.stack(gene_emb_dict[k]), dim=0), dim=0)
            .cpu()
            .numpy()
        )
    embs_df = pd.DataFrame(gene_emb_dict).T
    embs_df.index = [token_gene_dict[token] for token in embs_df.index]
    return embs_df


def plot_umap(embs_df, emb_dims, label, output_file, kwargs_dict, seed=0):
    only_embs_df = embs_df.iloc[:, :emb_dims]
    only_embs_df.index = pd.RangeIndex(0, only_embs_df.shape[0], name=None).astype(str)
    only_embs_df.columns = pd.RangeIndex(0, only_embs_df.shape[1], name=None).astype(
        str
    )
    vars_dict = {"embs": only_embs_df.columns}
    obs_dict = {"cell_id": list(only_embs_df.index), f"{label}": list(embs_df[label])}
    adata = anndata.AnnData(X=only_embs_df, obs=obs_dict, var=vars_dict)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sns.set(rc={"figure.figsize": (10, 10)}, font_scale=2.3)
    sns.set_style("white")
    default_kwargs_dict = {"size": 200}
    if kwargs_dict is not None:
        default_kwargs_dict.update(kwargs_dict)

    cats = set(embs_df[label])

    with plt.rc_context():
        ax = sc.pl.umap(adata, color=label, show=False, **default_kwargs_dict)
        ax.legend(
            markerscale=2,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            ncol=(1 if len(cats) <= 14 else 2 if len(cats) <= 30 else 3),
        )
        plt.show()
        plt.savefig(output_file, bbox_inches="tight")


def gen_heatmap_class_colors(labels, df):
    pal = sns.cubehelix_palette(
        len(Counter(labels).keys()),
        light=0.9,
        dark=0.1,
        hue=1,
        reverse=True,
        start=1,
        rot=-2,
    )
    lut = dict(zip(map(str, Counter(labels).keys()), pal))
    colors = pd.Series(labels, index=df.index).map(lut)
    return colors


def gen_heatmap_class_dict(classes, label_colors_series):
    class_color_dict_df = pd.DataFrame(
        {"classes": classes, "color": label_colors_series}
    )
    class_color_dict_df = class_color_dict_df.drop_duplicates(subset=["classes"])
    return dict(zip(class_color_dict_df["classes"], class_color_dict_df["color"]))


def make_colorbar(embs_df, label):
    labels = list(embs_df[label])

    cell_type_colors = gen_heatmap_class_colors(labels, embs_df)
    label_colors = pd.DataFrame(cell_type_colors, columns=[label])

    # create dictionary for colors and classes
    label_color_dict = gen_heatmap_class_dict(labels, label_colors[label])
    return label_colors, label_color_dict


def plot_heatmap(embs_df, emb_dims, label, output_file, kwargs_dict):
    sns.set_style("white")
    sns.set(font_scale=2)
    plt.figure(figsize=(15, 15), dpi=150)
    label_colors, label_color_dict = make_colorbar(embs_df, label)

    default_kwargs_dict = {
        "row_cluster": True,
        "col_cluster": True,
        "row_colors": label_colors,
        "standard_scale": 1,
        "linewidths": 0,
        "xticklabels": False,
        "yticklabels": False,
        "figsize": (15, 15),
        "center": 0,
        "cmap": "magma",
    }

    if kwargs_dict is not None:
        default_kwargs_dict.update(kwargs_dict)
    g = sns.clustermap(
        embs_df.iloc[:, 0:emb_dims].apply(pd.to_numeric), **default_kwargs_dict
    )

    plt.setp(g.ax_row_colors.get_xmajorticklabels(), rotation=45, ha="right")

    for label_color in list(label_color_dict.keys()):
        g.ax_col_dendrogram.bar(
            0, 0, color=label_color_dict[label_color], label=label_color, linewidth=0
        )

        g.ax_col_dendrogram.legend(
            title=f"{label}",
            loc="lower center",
            ncol=4,
            bbox_to_anchor=(0.5, 1),
            facecolor="white",
        )
    plt.show()
    logger.info(f"Output file: {output_file}")
    plt.savefig(output_file, bbox_inches="tight")


class EmbExtractor:
    valid_option_dict = {
        "model_type": {"Pretrained", "GeneClassifier", "CellClassifier"},
        "num_classes": {int},
        "emb_mode": {"cls", "cell", "gene"},
        "cell_emb_style": {"mean_pool"},
        "gene_emb_style": {"mean_pool"},
        "filter_data": {None, dict},
        "max_ncells": {None, int},
        "emb_layer": {-1, 0},
        "emb_label": {None, list},
        "labels_to_plot": {None, list},
        "forward_batch_size": {int},
        "token_dictionary_file": {None, str},
        "nproc": {int},
        "summary_stat": {None, "mean", "median", "exact_mean", "exact_median"},
        "genelist": {None, list},
    }

    def __init__(
        self,
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cls",
        cell_emb_style="mean_pool",
        gene_emb_style="mean_pool",
        filter_data=None,
        max_ncells=1000,
        emb_layer=-1,
        emb_label=None,
        labels_to_plot=None,
        forward_batch_size=100,
        nproc=4,
        summary_stat=None,
        token_dictionary_file=None,
        genelist=None,
    ):
        """
        Initialize embedding extractor.

        **[Par]ameters:**

        model_type : {"Pretrained", "GeneClassifier", "CellClassifier"}
            | Whether model is the pretrained Geneformer or a fine-tuned gene or cell classifier.
        num_classes : int
            | If model is a gene or cell classifier, specify number of classes it was trained to classify.
            | For the pretrained Geneformer model, number of classes is 0 as it is not a classifier.
        emb_mode : {"cls", "cell", "gene"}
            | Whether to output CLS, cell, or gene embeddings.
            | CLS embeddings are cell embeddings derived from the CLS token in the front of the rank value encoding.
        cell_emb_style : {"mean_pool"}
            | Method for summarizing cell embeddings if not using CLS token.
            | Currently only option is mean pooling of gene embeddings for given cell.
        gene_emb_style : "mean_pool"
            | Method for summarizing gene embeddings.
            | Currently only option is mean pooling of contextual gene embeddings for given gene.
        filter_data : None, dict
            | Default is to extract embeddings from all input data.
            | Otherwise, dictionary specifying .dataset column name and list of values to filter by.
        max_ncells : None, int
            | Maximum number of cells to extract embeddings from.
            | Default is 1000 cells randomly sampled from input data.
            | If None, will extract embeddings from all cells.
        emb_layer : {-1, 0}
            | Embedding layer to extract.
            | The last layer is most specifically weighted to optimize the given learning objective.
            | Generally, it is best to extract the 2nd to last layer to get a more general representation.
            | -1: 2nd to last layer
            | 0: last layer
        emb_label : None, list
            | List of column name(s) in .dataset to add as labels to embedding output.
        labels_to_plot : None, list
            | Cell labels to plot.
            | Shown as color bar in heatmap.
            | Shown as cell color in umap.
            | Plotting umap requires labels to plot.
        forward_batch_size : int
            | Batch size for forward pass.
        nproc : int
            | Number of CPU processes to use.
        summary_stat : {None, "mean", "median", "exact_mean", "exact_median"}
            | If exact_mean or exact_median, outputs only exact mean or median embedding of input data.
            | If mean or median, outputs only approximated mean or median embedding of input data.
            | Non-exact recommended if encountering memory constraints while generating goal embedding positions.
            | Non-exact is slower but more memory-efficient.
        token_dictionary_file : Path
            | Default is the Geneformer token dictionary
            | Path to pickle file containing token dictionary (Ensembl ID:token).

        **Examples:**

        .. code-block :: python

            >>> from geneformer import EmbExtractor
            >>> embex = EmbExtractor(model_type="CellClassifier",
            ...         num_classes=3,
            ...         emb_mode="cell",
            ...         filter_data={"cell_type":["cardiomyocyte"]},
            ...         max_ncells=1000,
            ...         emb_layer=-1,
            ...         emb_label=["disease", "cell_type"],
            ...         labels_to_plot=["disease", "cell_type"])

        """

        self.model_type = model_type
        self.num_classes = num_classes
        self.emb_mode = emb_mode
        self.cell_emb_style = cell_emb_style
        self.gene_emb_style = gene_emb_style
        self.filter_data = filter_data
        self.max_ncells = max_ncells
        self.emb_layer = emb_layer
        self.emb_label = emb_label
        self.labels_to_plot = labels_to_plot
        self.token_dictionary_file = token_dictionary_file
        self.forward_batch_size = forward_batch_size
        self.nproc = nproc
        self.genelist = genelist
        if (summary_stat is not None) and ("exact" in summary_stat):
            self.summary_stat = None
            self.exact_summary_stat = summary_stat
        else:
            self.summary_stat = summary_stat
            self.exact_summary_stat = None

        self.validate_options()

        with open(token_dictionary_file, "rb") as f:
            self.gene_token_dict = pickle.load(f)

        self.token_gene_dict = {v: k for k, v in self.gene_token_dict.items()}
        self.pad_token_id = self.gene_token_dict.get("<pad>")

    def validate_options(self):
        # confirm arguments are within valid options and compatible with each other
        for attr_name, valid_options in self.valid_option_dict.items():
            attr_value = self.__dict__[attr_name]
            if not isinstance(attr_value, (list, dict)):
                if attr_value in valid_options:
                    continue
            valid_type = False
            for option in valid_options:
                if (option in [int, list, dict, bool, str]) and isinstance(
                    attr_value, option
                ):
                    valid_type = True
                    break
            if valid_type:
                continue
            logger.error(
                f"Invalid option for {attr_name}. "
                f"Valid options for {attr_name}: {valid_options}"
            )
            raise

        if self.filter_data is not None:
            for key, value in self.filter_data.items():
                if not isinstance(value, list):
                    self.filter_data[key] = [value]
                    logger.warning(
                        "Values in filter_data dict must be lists. "
                        f"Changing {key} value to list ([{value}])."
                    )

    def extract_embs(
        self,
        model_directory,
        input_data_file,
        output_directory,
        output_prefix,
        output_torch_embs=False,
        cell_state=None,
        get_avg_attentions=False,
    ):
        """
        Extract embeddings from input data and save as results in output_directory.

        **[Par]ameters:**

        model_directory : Path
            | Path to directory containing model
        input_data_file : Path
            | Path to directory containing .dataset inputs
        output_directory : Path
            | Path to directory where embedding data will be saved as csv
        output_prefix : str
            | Prefix for output file
        output_torch_embs : bool
            | Whether or not to also output the embeddings as a tensor.
            | Note, if true, will output embeddings as both dataframe and tensor.
        cell_state : dict
            | Cell state key and value for state embedding extraction.

        **Examples:**

        .. code-block :: python

            >>> embs = embex.extract_embs("path/to/model",
            ...                           "path/to/input_data",
            ...                           "path/to/output_directory",
            ...                           "output_prefix")

        """

        filtered_input_data = pu.load_and_filter(
            self.filter_data, self.nproc, input_data_file
        )
        if cell_state is not None:
            filtered_input_data = pu.filter_by_dict(
                filtered_input_data, cell_state, self.nproc
            )
        downsampled_data = pu.downsample_and_sort(filtered_input_data, self.max_ncells)
        model = pu.load_model(
            self.model_type, self.num_classes, model_directory, mode="eval"
        )
        layer_to_quant = pu.quant_layers(model) + self.emb_layer
        embs, avg_attentions = get_embs(
            model=model,
            filtered_input_data=downsampled_data,
            emb_mode=self.emb_mode,
            layer_to_quant=layer_to_quant,
            pad_token_id=self.pad_token_id,
            forward_batch_size=self.forward_batch_size,
            token_gene_dict=self.token_gene_dict,
            summary_stat=self.summary_stat,
            get_avg_attentions=get_avg_attentions,
            genelist=self.genelist,
        )

        if self.emb_mode == "cell":
            if self.summary_stat is None:
                embs_df = label_cell_embs(embs, downsampled_data, self.emb_label)
            elif self.summary_stat is not None:
                embs_df = pd.DataFrame(embs.cpu().numpy()).T
        elif self.emb_mode == "gene":
            if self.summary_stat is None:
                embs_df = label_gene_embs(embs, downsampled_data, self.token_gene_dict)
            elif self.summary_stat is not None:
                embs_df = pd.DataFrame(embs).T
                embs_df.index = [self.token_gene_dict[token] for token in embs_df.index]
        elif self.emb_mode == "cls":
            embs_df = label_cell_embs(embs, downsampled_data, self.emb_label)

        # save embeddings to output_path
        if cell_state is None:
            output_path = (Path(output_directory) / output_prefix).with_suffix(".csv")
            embs_df.to_csv(output_path)

        if self.exact_summary_stat == "exact_mean":
            embs = embs.mean(dim=0)
            emb_dims = pu.get_model_emb_dims(model)
            embs_df = pd.DataFrame(
                embs_df[0 : emb_dims - 1].mean(axis="rows"),
                columns=[self.exact_summary_stat],
            ).T
        elif self.exact_summary_stat == "exact_median":
            embs = torch.median(embs, dim=0)[0]
            emb_dims = pu.get_model_emb_dims(model)
            embs_df = pd.DataFrame(
                embs_df[0 : emb_dims - 1].median(axis="rows"),
                columns=[self.exact_summary_stat],
            ).T

        if cell_state is not None:
            return embs, avg_attentions
        else:
            if output_torch_embs:
                return embs_df, embs, avg_attentions
            else:
                return embs_df, avg_attentions

    def get_state_embs(
        self,
        cell_states_to_model,
        model_directory,
        input_data_file,
        output_directory,
        output_prefix,
        output_torch_embs=True,
    ):
        """
        Extract exact mean or exact median cell state embedding positions from input data and save as results in output_directory.

        **[Par]ameters:**

        cell_states_to_model : None, dict
            | Cell states to model if testing perturbations that achieve goal state change.
            | Four-item dictionary with keys: state_key, start_state, goal_state, and alt_states
            | state_key: key specifying name of column in .dataset that defines the start/goal states
            | start_state: value in the state_key column that specifies the start state
            | goal_state: value in the state_key column taht specifies the goal end state
            | alt_states: list of values in the state_key column that specify the alternate end states
            | For example:
            |      {"state_key": "disease",
            |      "start_state": "dcm",
            |      "goal_state": "nf",
            |      "alt_states": ["hcm", "other1", "other2"]}
        model_directory : Path
            | Path to directory containing model
        input_data_file : Path
            | Path to directory containing .dataset inputs
        output_directory : Path
            | Path to directory where embedding data will be saved as csv
        output_prefix : str
            | Prefix for output file
        output_torch_embs : bool
            | Whether or not to also output the embeddings as a tensor.
            | Note, if true, will output embeddings as both dataframe and tensor.

        **Outputs**

        | Outputs state_embs_dict for use with in silico perturber.
        | Format is dictionary of embedding positions of each cell state to model shifts from/towards.
        | Keys specify each possible cell state to model.
        | Values are target embedding positions as torch.tensor.
        | For example:
        |      {"nf": emb_nf,
        |      "hcm": emb_hcm,
        |      "dcm": emb_dcm,
        |      "other1": emb_other1,
        |      "other2": emb_other2}
        """

        pu.validate_cell_states_to_model(cell_states_to_model)
        valid_summary_stats = ["exact_mean", "exact_median"]
        if self.exact_summary_stat not in valid_summary_stats:
            logger.error(
                "For extracting state embs, summary_stat in EmbExtractor "
                f"must be set to option in {valid_summary_stats}"
            )
            raise

        if self.emb_label is not None:
            logger.error(
                "For extracting state embs, emb_label should be None since labels are based on state embs dict keys."
            )
            raise

        state_embs_dict = dict()
        state_key = cell_states_to_model["state_key"]
        for k, v in cell_states_to_model.items():
            if k == "state_key":
                continue
            elif (k == "start_state") or (k == "goal_state"):
                state_embs_dict[v] = self.extract_embs(
                    model_directory,
                    input_data_file,
                    output_directory,
                    output_prefix,
                    output_torch_embs,
                    cell_state={state_key: v},
                )
            else:  # k == "alt_states"
                for alt_state in v:
                    state_embs_dict[alt_state] = self.extract_embs(
                        model_directory,
                        input_data_file,
                        output_directory,
                        output_prefix,
                        output_torch_embs,
                        cell_state={state_key: alt_state},
                    )

        output_path = (Path(output_directory) / output_prefix).with_suffix(".pkl")
        with open(output_path, "wb") as fp:
            pickle.dump(state_embs_dict, fp)

        return state_embs_dict

    def plot_embs(
        self,
        embs,
        plot_style,
        output_directory,
        output_prefix,
        max_ncells_to_plot=1000,
        kwargs_dict=None,
    ):
        """
        Plot embeddings, coloring by provided labels.

        **Parameters:**

        embs : pandas.core.frame.DataFrame
            | Pandas dataframe containing embeddings output from extract_embs
        plot_style : str
            | Style of plot: "heatmap" or "umap"
        output_directory : Path
            | Path to directory where plots will be saved as pdf
        output_prefix : str
            | Prefix for output file
        max_ncells_to_plot : None, int
            | Maximum number of cells to plot.
            | Default is 1000 cells randomly sampled from embeddings.
            | If None, will plot embeddings from all cells.
        kwargs_dict : dict
            | Dictionary of kwargs to pass to plotting function.

        **Examples:**

        .. code-block :: python

            >>> embex.plot_embs(embs=embs,
            ...                 plot_style="heatmap",
            ...                 output_directory="path/to/output_directory",
            ...                 output_prefix="output_prefix")

        """

        if plot_style not in ["heatmap", "umap"]:
            logger.error(
                "Invalid option for 'plot_style'. Valid options: {'heatmap','umap'}"
            )
            raise

        if (plot_style == "umap") and (self.labels_to_plot is None):
            logger.error("Plotting UMAP requires 'labels_to_plot'. ")
            raise

        if max_ncells_to_plot is not None:
            if max_ncells_to_plot > self.max_ncells:
                max_ncells_to_plot = self.max_ncells
                logger.warning(
                    "max_ncells_to_plot must be <= max_ncells. "
                    f"Changing max_ncells_to_plot to {self.max_ncells}."
                )
            elif max_ncells_to_plot < self.max_ncells:
                embs = embs.sample(max_ncells_to_plot, axis=0)

        if self.emb_label is None:
            label_len = 0
        else:
            label_len = len(self.emb_label)

        emb_dims = embs.shape[1] - label_len

        if self.emb_label is None:
            emb_labels = None
        else:
            emb_labels = embs.columns[emb_dims:]

        if plot_style == "umap":
            for label in self.labels_to_plot:
                if label not in emb_labels:
                    logger.warning(
                        f"Label {label} from labels_to_plot "
                        f"not present in provided embeddings dataframe."
                    )
                    continue
                output_prefix_label = output_prefix + f"_umap_{label}"
                output_file = (
                    Path(output_directory) / output_prefix_label
                ).with_suffix(".pdf")
                plot_umap(embs, emb_dims, label, output_file, kwargs_dict)

        if plot_style == "heatmap":
            for label in self.labels_to_plot:
                if label not in emb_labels:
                    logger.warning(
                        f"Label {label} from labels_to_plot "
                        f"not present in provided embeddings dataframe."
                    )
                    continue
                output_prefix_label = output_prefix + f"_heatmap_{label}"
                output_file = (
                    Path(output_directory) / output_prefix_label
                ).with_suffix(".pdf")
                plot_heatmap(embs, emb_dims, label, output_file, kwargs_dict)


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


def parse_args(par):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rna", type=str, help="Path to the input RNA data in h5ad format."
    )
    parser.add_argument(
        "--rna_all", type=str, help="Path to the input RNA all data in h5ad format."
    )
    parser.add_argument(
        "--atac", type=str, help="Path to the input ATAC data in h5ad format."
    )
    parser.add_argument(
        "--prediction", type=str, help="Path to the output prediction in h5ad format."
    )
    parser.add_argument(
        "--score", type=str, help="Path to the output score in h5ad format."
    )
    parser.add_argument("--layer", type=str)
    parser.add_argument("--temp_dir", type=str)
    parser.add_argument("--tf_all", type=str)
    parser.add_argument("--skeleton", type=str)
    parser.add_argument("--apply_skeleton", action="store_true")
    parser.add_argument("--apply_tf", action="store_true")
    parser.add_argument("--max_n_links", type=int)
    parser.add_argument("--reg_type", type=str)
    parser.add_argument("--num_workers", type=int)
    parser.add_argument("--regulators_consensus", type=str)
    parser.add_argument("--evaluation_data", type=str)
    parser.add_argument("--ws_consensus", type=str)
    parser.add_argument("--ws_distance_background", type=str)

    args = parser.parse_args()
    for k, v in vars(args).items():
        if v is not None:
            par[k] = v
    return par


par = parse_args(par)

## VIASH END

n_processors = os.cpu_count()

## GETTING DATA

adata = sc.read(par["rna"], backed="r")

if adata.uns["dataset_id"] in ["replogle", "xaira_HCT116", "xaira_HEK293T"]:
    train_perturbs = adata.obs["perturbation"].unique()
    tf_all = np.loadtxt(par["tf_all"], dtype=str)
    train_perturbs = np.intersect1d(tf_all, train_perturbs)
    train_perturbs = train_perturbs[:100]  # limit to 100 perturbations
    mask = adata.obs["perturbation"].isin(train_perturbs)
    adata = adata[mask].to_memory()
elif adata.uns["dataset_id"] in ["parsebioscience"]:
    train_perturbs = adata.obs["perturbation"].unique()
    train_perturbs = train_perturbs[:10]
    mask = adata.obs["perturbation"].isin(train_perturbs)
    adata = adata[mask].to_memory()
else:
    adata = adata.to_memory()

if adata.raw is not None and adata.raw.X.shape[1] != adata.X.shape[1]:
    print("removing raw")
    del adata.raw
if adata.layers is not None and "counts" in adata.layers:
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]

if adata[0].X.sum() != int(adata[0].X.sum()):
    print("WARNING: you are not using count data")
    print("reverting logp1")
    adata.X = csr_matrix(np.power(adata.X.todense(), 2) - 1)

adata.var["symbol"] = adata.var.index
adata.var["ensembl_id"] = adata.var["gene_ids"].values
dataset_id = adata.uns["dataset_id"] if "dataset_id" in adata.uns else par["dataset_id"]


### LOADING MODEL AND PREDICTION

print(f">>> Getting settings for model '{par['model']}'...", flush=True)
model_split = par["model"].split("-")
model_details = {
    "dataset": "104M",
    "input_size": 768,
}

print(">>> Getting model dictionary files...", flush=True)
gene_mapping_file = hf_hub_download(
    repo_id="ctheodoris/Geneformer",
    subfolder="geneformer",
    filename=f"ensembl_mapping_dict_gc{model_details['dataset']}.pkl",
)
gene_median = hf_hub_download(
    repo_id="ctheodoris/Geneformer",
    subfolder="geneformer",
    filename=f"gene_median_dictionary_gc{model_details['dataset']}.pkl",
)
gene_name_id = hf_hub_download(
    repo_id="ctheodoris/Geneformer",
    subfolder="geneformer",
    filename=f"gene_name_id_dict_gc{model_details['dataset']}.pkl",
)
token = hf_hub_download(
    repo_id="ctheodoris/Geneformer",
    subfolder="geneformer",
    filename=f"token_dictionary_gc{model_details['dataset']}.pkl",
)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

#### PREDICTION

if "cell_type" not in adata.obs:
    adata.obs["cell_type"] = "dummy_cell_type"
    par["how"] = "most var within"


print(">>> Creating working directory...", flush=True)
work_dir = TemporaryDirectory()
input_dir = os.path.join(work_dir.name, "input")
os.makedirs(input_dir)
tokenized_dir = os.path.join(work_dir.name, "tokenized")
os.makedirs(tokenized_dir)
embedding_dir = os.path.join(work_dir.name, "embedding")
os.makedirs(embedding_dir)
print(f"Working directory: '{work_dir.name}'", flush=True)

print(">>> Preparing data...", flush=True)
adata.var["ensembl_id"] = adata.var_names
adata.obs["n_counts"] = np.ravel(adata.X.sum(axis=1))
adata.write_h5ad(os.path.join(input_dir, "input.h5ad"))
print(adata)


# Function to try parallel execution and fall batch to a single processor if it fails
def tryParallelFunction(fun, label):
    try:
        fun(nproc=n_processors)
    except RuntimeError as e:
        # Retry with nproc=1 if error message contains "One of the subprocesses has abruptly died"
        if "subprocess" in str(e) and "died" in str(e):
            print(f"{label} failed. Error message: {e}", flush=True)
            print("Retrying with nproc=1", flush=True)
            fun(nproc=1)
        else:
            raise e


print(">>> Tokenizing data...", flush=True)


def tokenize_data(nproc):
    tokenizer = TranscriptomeTokenizer(
        nproc=nproc,
        model_input_size=model_details["input_size"],
        gene_median_file=gene_median,
        token_dictionary_file=token,
        gene_mapping_file=gene_mapping_file,
    )

    tokenizer.tokenize_data(
        "/tmp/geneformer/", tokenized_dir, "tokenized", file_format="h5ad"
    )


model_files = {
    "model": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=par["model"],
        filename="model.safetensors",
    ),
    "config": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=par["model"],
        filename="config.json",
    ),
}

model_dir = os.path.dirname(model_files["model"])

with open(gene_mapping_file, "rb") as f:
    gene_mapping_dict = pickle.load(f)


def compute_geneformer_network(
    adata,
    forward_batch_size=4,
    max_ncells=1000,
):
    adata = adata[
        :,
        adata.var.index.isin(
            [u for u in gene_mapping_dict.keys() if u in adata.var.index]
        ),
    ]
    adata.obs["n_counts"] = adata.X.sum(1)
    # Create the geneformer folder if it doesn't exist
    geneformer_folder = "/tmp/geneformer"
    if not os.path.exists(geneformer_folder):
        os.makedirs(geneformer_folder)
    adata.write_h5ad("/tmp/geneformer/to_token.h5ad")

    genelist = [gene_mapping_dict[u] for u in adata.var.index]

    tokenized_data_path = "/tmp/geneformer/tokenized_data.dataset"
    if os.path.exists(tokenized_data_path):
        shutil.rmtree(tokenized_data_path)

    tryParallelFunction(tokenize_data, "Tokenizing data")

    embex = EmbExtractor(
        model_type="Pretrained",  # CellClassifier
        max_ncells=max_ncells,
        emb_mode="cell",
        forward_batch_size=forward_batch_size,
        nproc=n_processors,
        genelist=genelist,
        token_dictionary_file=token,
    )
    embs, avg_attentions = embex.extract_embs(
        model_dir,  # ../../../Geneformer/fine_tuned_models/gf-6L-30M-i2048_CellClassifier_cardiomyopathies_220224/
        os.path.join(tokenized_dir, "tokenized.dataset"),
        embedding_dir,
        "embedding",
        get_avg_attentions=True,
    )
    return adata.copy(), avg_attentions.to("cpu").numpy()


for i, cell_type in enumerate(adata.obs["cell_type"].unique()):
    print(cell_type)
    subadata = adata[adata.obs["cell_type"] == cell_type].copy()
    subadata = subadata[
        : par["max_cells"], subadata.X.sum(0).A1.argsort()[::-1][: par["num_genes"]]
    ]
    subadata, net = compute_geneformer_network(
        subadata,
        forward_batch_size=par["batch_size"],
    )
    gene_names = subadata.var["symbol"].values
    print(net.shape, net.sum(), len(gene_names), gene_names[:10])
    net = efficient_melting(net.T, gene_names, symmetric=False)
    net = net[net["weight"] != 0]

    # - subset to TFs
    tfs = set(np.loadtxt(par["tf_all"], dtype=str))

    net = net[net["source"].isin(tfs)]

    net = net.sort_values(
        by="weight", ascending=False, key=abs
    )[
        : 2 * par["max_n_links"]
    ]  # i set this to double of allowed link just in case the symmetry exists. metrics will take care of this
    net["cell_type"] = cell_type
    if i == 0:
        net_all = net
    else:
        net_all = pd.concat([net_all, net], axis=0)


print(f"Writing results to {par['prediction']}")
net_all = net_all.astype(str)
output = ad.AnnData(
    X=None,
    uns={
        "method_id": "scprint",
        "dataset_id": dataset_id,
        "prediction": net_all[["source", "target", "weight", "cell_type"]],
    },
)
output.write(par["prediction"])
output.write(par["prediction"])
