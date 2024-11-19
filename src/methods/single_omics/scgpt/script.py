import copy
import json
import os
from pathlib import Path
import sys
import warnings
import subprocess

import torch
from anndata import AnnData
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import tqdm
import os
# import gseapy as gp
# from gears import PertData, GEARS

from scipy.sparse import issparse
import scipy as sp
import numpy as np
from einops import rearrange
from torch.nn.functional import softmax
from tqdm import tqdm

from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)
import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils import set_seed 
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.preprocess import Preprocessor

os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')

# torch._dynamo.config.optimize_ddp = False

## VIASH START
par = {
  'multiomics_rna': 'resources_test/grn-benchmark/multiomics_rna.h5ad',
  'tf_all': 'resources_test/prior/tf_all.csv',
  'prediction': 'output/prediction_scgpt.csv',
  'max_n_links': 50000,
  'model_file': 'resources_test/supplementary/finetuned_scGPT_adamson/best_model.pt',
  'model_config_file': 'resources_test/supplementary/finetuned_scGPT_adamson/args.json',
  'vocab_file': 'resources_test/supplementary/finetuned_scGPT_adamson/vocab.json',
  'n_bins': 51,
  'batch_size': 16,
  'condition': 'cell_type'
}
## VIASH END

# Download datasets 
par['model_file'] = f"{par['temp_dir']}/best_model.pt"
par['model_config_file'] = f"{par['temp_dir']}/args.json"
par['vocab_file'] = f"{par['temp_dir']}/vocab.json"


import requests
def download_file(output_file, url):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
        print(f"File downloaded successfully and saved to {output_file}")
    else:
        print(f"Failed to download file. HTTP status code: {response.status_code}")
download_file(par['model_file'], 'https://drive.google.com/uc?export=download&id=1CPVtpWUJ2nkI9jGignlHLcefBe6Gk-F9')
download_file(par['vocab_file'], 'https://drive.google.com/file/d/1Qzb6Y9UB342a2QxmY-BCubSvcmYZ5jw3/view?usp=drive_link')
download_file(par['model_config_file'], 'https://drive.google.com/file/d/1VwPGHuSorVAXyTreMFI1yzMougtUDeUt/view?usp=drive_link')


# os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:50"
initial_memory = torch.cuda.memory_allocated()
def monitor_memory():
  used_memory = torch.cuda.memory_allocated()
  data_moved = used_memory - initial_memory
  print(f"Data moved to GPU: {data_moved} bytes")

# Load list of putative TFs
tf_all = np.loadtxt(par['tf_all'], dtype=str)

set_seed(42)
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]

mask_value = -1
pad_value = -2
batch_size = par['batch_size']
num_attn_layers = 11 
n_input_bins = par['n_bins']


model_config_file = par['model_config_file']
model_file = par['model_file']
vocab_file = par['vocab_file']

vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)
print('Reading the pretrained model model')


# Retrieve model parameters from config files
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
print(
    f"Resume model from {model_file}, the model args will override the "
    f"config {model_config_file}."
)
embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]

gene2idx = vocab.get_stoi()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

ntokens = len(vocab)  # size of vocabulary
model = TransformerModel(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    vocab=vocab,
    pad_value=pad_value,
    n_input_bins=n_input_bins,
    use_fast_transformer=True,
    domain_spec_batchnorm = "batchnorm"
)

try:
    model.load_state_dict(torch.load(model_file))
    print(f"Loading all model params from {model_file}")
except:
    # only load params that are in the model and match the size
    model_dict = model.state_dict()
    pretrained_dict = torch.load(model_file)
    pretrained_dict = {
        k: v
        for k, v in pretrained_dict.items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    for k, v in pretrained_dict.items():
        print(f"Loading params {k} with shape {v.shape}")
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict)

model.to(device)
monitor_memory()


print('Process rna-seq file')
import scanpy as sc 
adata = sc.read(par['multiomics_rna'])
adata.X = adata.X.todense()
adata.obs["celltype"] = adata.obs["cell_type"].astype("category")
adata.obs["str_batch"] = adata.obs["donor_id"].astype(str)
data_is_raw = False

adata.var["id_in_vocab"] = [1 if gene in vocab else -1 for gene in adata.var.index]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=3,  # step 1
    filter_cell_by_counts=False,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=data_is_raw,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg= False,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
    binning=n_input_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)
preprocessor(adata, batch_key="str_batch")

# print('Subsetting to HVGs')
# sc.pp.highly_variable_genes(
#     adata,
#     layer=None,
#     n_top_genes=n_hvg,
#     batch_key="str_batch",
#     flavor="seurat_v3" if data_is_raw else "cell_ranger",
#     subset=False,
# )
# adata = adata[:, adata.var["highly_variable"]].copy()


input_layer_key = "X_binned"
all_counts = (
    adata.layers[input_layer_key].A
    if issparse(adata.layers[input_layer_key])
    else adata.layers[input_layer_key]
)
genes = adata.var.index.tolist()
gene_ids = np.array(vocab(genes), dtype=int)

batch_size = batch_size
tokenized_all = tokenize_and_pad_batch(
    all_counts,
    gene_ids,
    max_len=len(genes)+1,
    vocab=vocab,
    pad_token=pad_token,
    pad_value=pad_value,
    append_cls=True,  # append <cls> token at the beginning
    include_zero_gene=True,
)


all_gene_ids, all_values = tokenized_all["genes"], tokenized_all["values"]


src_key_padding_mask = all_gene_ids.eq(vocab[pad_token])

condition_ids = np.array(adata.obs[par['condition']].tolist())

torch.cuda.empty_cache()
dict_sum_condition = {}
print('Extract gene gene links from attention layer')
model.eval()
monitor_memory()
with torch.no_grad(), torch.cuda.amp.autocast(enabled=True):
    M = all_gene_ids.size(1)
    N = all_gene_ids.size(0)
    device = next(model.parameters()).device
    for i in tqdm(range(0, N, batch_size)):
        batch_size = all_gene_ids[i : i + batch_size].size(0)
        outputs = np.zeros((batch_size, M, M), dtype=np.float32)
        # Replicate the operations in model forward pass
        src_embs = model.encoder(torch.tensor(all_gene_ids[i : i + batch_size], dtype=torch.long).to(device))
        # monitor_memory()
        val_embs = model.value_encoder(torch.tensor(all_values[i : i + batch_size], dtype=torch.float).to(device))
        total_embs = src_embs + val_embs
        total_embs = model.bn(total_embs.permute(0, 2, 1)).permute(0, 2, 1)
        # Send total_embs to attention layers for attention operations
        # Retrieve the output from second to last layer
        for layer in model.transformer_encoder.layers[:num_attn_layers]:
            total_embs = layer(total_embs, src_key_padding_mask=src_key_padding_mask[i : i + batch_size].to(device))
        # Send total_embs to the last layer in flash-attn
        # https://github.com/HazyResearch/flash-attention/blob/1b18f1b7a133c20904c096b8b222a0916e1b3d37/flash_attn/flash_attention.py#L90
        qkv = model.transformer_encoder.layers[num_attn_layers].self_attn.Wqkv(total_embs)
        # Retrieve q, k, and v from flast-attn wrapper
        qkv = rearrange(qkv, 'b s (three h d) -> b s three h d', three=3, h=8)
        q = qkv[:, :, 0, :, :]
        k = qkv[:, :, 1, :, :]
        v = qkv[:, :, 2, :, :]
        # https://towardsdatascience.com/illustrated-self-attention-2d627e33b20a
        # q = [batch, gene, n_heads, n_hid]
        # k = [batch, gene, n_heads, n_hid]
        # attn_scores = [batch, n_heads, gene, gene]
        attn_scores = q.permute(0, 2, 1, 3) @ k.permute(0, 2, 3, 1)
        # Rank normalization by row
        attn_scores = attn_scores.reshape((-1, M))
        order = torch.argsort(attn_scores, dim=1)
        rank = torch.argsort(order, dim=1)
        attn_scores = rank.reshape((-1, 8, M, M))/M
        # Rank normalization by column
        attn_scores = attn_scores.permute(0, 1, 3, 2).reshape((-1, M))
        order = torch.argsort(attn_scores, dim=1)
        rank = torch.argsort(order, dim=1)
        attn_scores = (rank.reshape((-1, 8, M, M))/M).permute(0, 1, 3, 2)
        # Average 8 attention heads
        attn_scores = attn_scores.mean(1)
        outputs = attn_scores.detach().cpu().numpy()
        for index in range(batch_size):
            # Keep track of sum per condition
            c = condition_ids[i : i + batch_size][index]
            if c not in dict_sum_condition:
                dict_sum_condition[c] = np.zeros((M, M), dtype=np.float32)
            else:
                dict_sum_condition[c] += outputs[index, :, :]
print('Average across groups of cell types')
groups = adata.obs.groupby([par['condition']]).groups
dict_sum_condition_mean = dict_sum_condition.copy()
for i in groups.keys():
    dict_sum_condition_mean[i] = dict_sum_condition_mean[i]/len(groups[i])
mean_grn = np.array(list(dict_sum_condition_mean.values())).mean(axis=0)
print('Subset only for TFs')
gene_vocab_idx = all_gene_ids[0].clone().detach().cpu().numpy()
gene_names = vocab.lookup_tokens(gene_vocab_idx)

print('Format as df, melt, and subset')
net = pd.DataFrame(mean_grn, columns=gene_names, index=gene_names)
net = net.iloc[1:, 1:]

tf_all = np.intersect1d(tf_all, gene_names)
net = net[tf_all]

net_melted = net.reset_index()  # Move index to a column for melting
net_melted = pd.melt(net_melted, id_vars=net_melted.columns[0], var_name='target', value_name='weight')
net_melted.rename(columns={net_melted.columns[0]: 'source'}, inplace=True)

net_melted.to_csv(par['prediction'])