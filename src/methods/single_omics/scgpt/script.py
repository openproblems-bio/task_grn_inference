import copy
import json
import os
import sys
import warnings
import torch
import scanpy as sc
import numpy as np
import pandas as pd
import tqdm

from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed

os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')

## VIASH START
par = {
  'rna': 'resources/grn_benchmark/inference_data/replogle_rna.h5ad',
  "tf_all": 'resources/grn_benchmark/prior/tf_all.csv',
  'prediction': 'output/scgpt_test.h5ad',
  'temp_dir': 'output/scgpt',
  "model": "scGPT_human",
  'num_workers': 10
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
      'util_dir': 'src/utils',
      'helper_dir': 'src/methods/single_omics/scgpt',
    }
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])
from util import process_links
# from helper import format_data


if par["model"] is None:
    print(f"\n>>> Downloading '{par['model_name']}' model...", flush=True)
    model_drive_ids = {
        "scGPT_human": "1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
        "scGPT_CP": "1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB",
    }
    drive_path = (
        f"https://drive.google.com/drive/folders/{model_drive_ids[par['model_name']]}"
    )
    model_temp = tempfile.TemporaryDirectory()
    model_dir = model_temp.name
    print(f"Downloading from '{drive_path}'", flush=True)
    gdown.download_folder(drive_path, output=model_dir, quiet=True)
else:
    if os.path.isdir(par["model"]):
        print(f"\n>>> Using model directory...", flush=True)
        model_temp = None
        model_dir = par["model"]
    else:
        model_temp = tempfile.TemporaryDirectory()
        model_dir = model_temp.name

        if zipfile.is_zipfile(par["model"]):
            print(f"\n>>> Extracting model from .zip...", flush=True)
            print(f".zip path: '{par['model']}'", flush=True)
            with zipfile.ZipFile(par["model"], "r") as zip_file:
                zip_file.extractall(model_dir)
        elif tarfile.is_tarfile(par["model"]) and par["model"].endswith(
            ".tar.gz"
        ):
            print(f"\n>>> Extracting model from .tar.gz...", flush=True)
            print(f".tar.gz path: '{par['model']}'", flush=True)
            with tarfile.open(par["model"], "r:gz") as tar_file:
                tar_file.extractall(model_dir)
                model_dir = os.path.join(model_dir, os.listdir(model_dir)[0])
        else:
            raise ValueError(
                f"The 'model' argument should be a directory a .zip file or a .tar.gz file"
            )

print(f"Model directory: '{model_dir}'", flush=True)


set_seed(42)
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
n_bins = 51
mask_value = -1
pad_value = -2
n_input_bins = n_bins

model_config_file = model_dir / "args.json"
model_file = model_dir / "best_model.pt"
vocab_file = model_dir / "vocab.json"

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