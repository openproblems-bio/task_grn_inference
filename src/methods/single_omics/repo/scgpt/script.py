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

import requests
import gdown


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
  'rna': 'resources_test/grn-benchmark/multiomics_rna.h5ad',
  'tf_all': 'resources_test/grn_benchmark/prior/tf_all.csv',
  'prediction': 'output/prediction_scgpt.csv',
  'max_n_links': 50000,
  'model_file': 'resources_test/supplementary/finetuned_scGPT_adamson/best_model.pt',
  'model_config_file': 'resources_test/supplementary/finetuned_scGPT_adamson/args.json',
  'vocab_file': 'resources_test/supplementary/finetuned_scGPT_adamson/vocab.json',
  'n_bins': 51,
  'batch_size': 16,
  'condition': 'cell_type',
  'temp_dir': 'output'
}
## VIASH END
os.makedirs(par['temp_dir'], exist_ok=True)
# Download datasets 
par['model_file'] = f"{par['temp_dir']}/best_model.pt"
par['model_config_file'] = f"{par['temp_dir']}/args.json"
par['vocab_file'] = f"{par['temp_dir']}/vocab.json"

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
download_file(par['vocab_file'], 'https://docs.google.com/uc?export=download&id=1Qzb6Y9UB342a2QxmY-BCubSvcmYZ5jw3')
download_file(par['model_config_file'], 'https://docs.google.com/uc?export=download&id=1VwPGHuSorVAXyTreMFI1yzMougtUDeUt')


gdown.download("https://drive.google.com/uc?id=1CPVtpWUJ2nkI9jGignlHLcefBe6Gk-F9", par['model_file'], quiet=False)

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
adata = sc.read(par['rna'])
adata.X = adata.X.todense()
adata.obs["celltype"] = adata.obs["cell_type"].astype("str")
adata.obs["str_batch"] = adata.obs["donor_id"].astype(str)
data_is_raw = True
preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=3,  # step 1
    filter_cell_by_counts=False,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=data_is_raw,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg=False,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
    binning=51,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)
preprocessor(adata, batch_key="batch")

# Retrieve the data-independent gene embeddings from scGPT
gene_ids = np.array([id for id in gene2idx.values()])
gene_embeddings = model.encoder(torch.tensor(gene_ids, dtype=torch.long).to(device))
gene_embeddings = gene_embeddings.detach().cpu().numpy()

# Filter on the intersection between the Immune Human HVGs found in step 1.2 and scGPT's 30+K foundation model vocab
gene_embeddings = {gene: gene_embeddings[i] for i, gene in enumerate(gene2idx.keys()) if gene in adata.var.index.tolist()}
print('Retrieved gene embeddings for {} genes.'.format(len(gene_embeddings)))


embed = GeneEmbedding(gene_embeddings)


# Perform Louvain clustering with desired resolution; here we specify resolution=40
gdata = embed.get_adata(resolution=40)
# Retrieve the gene clusters
metagenes = embed.get_metagenes(gdata)

# Obtain the set of gene programs from clusters with #genes >= 5
mgs = dict()
for mg, genes in metagenes.items():
    if len(genes) > 4:
        mgs[mg] = genes

        

net = net_melted
net['weight'] = net['weight'].astype(str)
output = ad.AnnData(X=None, uns={"method_id": par['method_id'], "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
output.write(par['prediction'])