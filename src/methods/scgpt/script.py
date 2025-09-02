import copy
import json
import os
from pathlib import Path
import sys
import warnings

import torch
from anndata import AnnData
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import tqdm
import gdown


sys.path.insert(0, "../")
import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed

os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')


set_seed(42)
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
n_hvg = 1200
n_bins = 51
mask_value = -1
pad_value = -2
n_input_bins = n_bins

## VIASH START
par = {
  'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
  'tf_all': 'resources_test/prior/tf_all.csv',
  'prediction': 'output/scgpt/prediction.h5ad',
  'max_n_links': 50000,
#   'n_bins': 51,
  'batch_size': 16,
#   'condition': 'cell_type',
  'temp_dir': 'output/scgpt/'
}
## VIASH END
os.makedirs(par['temp_dir'], exist_ok=True)
# Download datasets 
model_file = f"{par['temp_dir']}/best_model.pt"
model_config_file = f"{par['temp_dir']}/args.json"
vocab_file = f"{par['temp_dir']}/vocab.json"

if os.path.exists(vocab_file) is False:
    def download_file(output_file, url):
        import requests
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(output_file, "wb") as f:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
            print(f"File downloaded successfully and saved to {output_file}")
        else:
            print(f"Failed to download file. HTTP status code: {response.status_code}")
    download_file(vocab_file, 'https://docs.google.com/uc?export=download&id=1Qzb6Y9UB342a2QxmY-BCubSvcmYZ5jw3')
if os.path.exists(model_config_file) is False:
    download_file(model_config_file, 'https://docs.google.com/uc?export=download&id=1VwPGHuSorVAXyTreMFI1yzMougtUDeUt')
if os.path.exists(model_file) is False:
    gdown.download("https://drive.google.com/uc?id=1CPVtpWUJ2nkI9jGignlHLcefBe6Gk-F9", model_file, quiet=False)


vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

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
)

try:
    model.load_state_dict(torch.load(model_file))
    # model.load_state_dict(torch.load(model_file, map_location=torch.device('cpu')))
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



