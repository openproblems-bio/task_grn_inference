import copy
import gc
import json
import time
from pathlib import Path

import numpy as np
import torch
import wandb
from einops import rearrange
from scgpt import logger
from scgpt.loss import masked_mse_loss
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.trainer import evaluate as scgpt_validate
from scgpt.trainer import prepare_data, prepare_dataloader
from scgpt.trainer import train as scgpt_train
from scgpt.utils import add_file_handler, set_seed
from scipy.sparse import issparse
from sklearn.model_selection import train_test_split
from torch import nn
from torch.nn import functional as F
from torch.nn.functional import softmax
from tqdm import tqdm


def prepare_model(
    model_dir="../save/scGPT_human",
    special_tokens=["<pad>", "<cls>", "<eoc>"],
    n_bins=51,
    pad_value=-2,
):
    # Specify model path; here we load the scGPT blood model fine-tuned on adamson
    model_dir = Path(model_dir)
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
        n_input_bins=n_bins,
        use_fast_transformer=True,
    )

    try:
        model.load_state_dict(torch.load(model_file))
        # print(f"Loading all model params from {model_file}")
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
            # print(f"Loading params {k} with shape {v.shape}")
            model_dict.update(pretrained_dict)
            model.load_state_dict(model_dict)

    model.to(device)
    return model, vocab


def prepare_dataset(
    subadata, vocab, data_is_raw=True, pad_token="<pad>", n_bins=51, pad_value=-2
):
    subadata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in subadata.var.index.astype(str)
    ]
    subadata = subadata[:, subadata.var["id_in_vocab"] >= 0]

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
        binning=n_bins,  # 6. whether to bin the raw data and to what number of bins
        result_binned_key="X_binned",  # the key in adata.layers to store the binned data
    )
    preprocessor(subadata, batch_key="str_batch")

    input_layer_key = "X_binned"
    all_counts = (
        subadata.layers[input_layer_key].A
        if issparse(subadata.layers[input_layer_key])
        else subadata.layers[input_layer_key]
    )
    genes = subadata.var.index.tolist()
    gene_ids = np.array(vocab(genes), dtype=int)
    tokenized_all = tokenize_and_pad_batch(
        all_counts,
        gene_ids,
        max_len=len(genes) + 1,
        vocab=vocab,
        pad_token=pad_token,
        pad_value=pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=True,
    )

    all_gene_ids, all_values = tokenized_all["genes"], tokenized_all["values"]
    src_key_padding_mask = all_gene_ids.eq(vocab[pad_token])
    return all_gene_ids, all_values, src_key_padding_mask, subadata


def generate_embedding(model, vocab, adata, batch_size=10):
    torch.cuda.empty_cache()
    all_gene_ids, all_values, src_key_padding_mask, adata = prepare_dataset(
        adata, vocab
    )
    model.eval()
    with torch.no_grad(), torch.cuda.amp.autocast(enabled=True):
        cell_embeddings = model.encode_batch(
            all_gene_ids,
            all_values,
            src_key_padding_mask,
            batch_size=batch_size,
            time_step=0,
        )
        cell_embeddings = cell_embeddings / np.linalg.norm(
            cell_embeddings, axis=1, keepdims=True
        )
    adata.obsm["scgpt_emb"] = cell_embeddings
    return adata


def generate_grn(model, vocab, adata, batch_size=10, num_attn_layers=11):
    torch.cuda.empty_cache()
    dict_sum_condition = {}
    all_gene_ids, all_values, src_key_padding_mask, subadata = prepare_dataset(
        adata, vocab
    )
    # Use this argument to specify which layer to extract the attention weights from
    # Default to 11, extraction from the last (12th) layer. Note that index starts from 0
    model.eval()
    with torch.no_grad(), torch.cuda.amp.autocast(enabled=True):
        M = all_gene_ids.size(1)
        N = all_gene_ids.size(0)
        device = next(model.parameters()).device
        for i in tqdm(range(0, N, batch_size)):
            batch_size = all_gene_ids[i : i + batch_size].size(0)
            outputs = np.zeros((batch_size, M, M), dtype=np.float32)
            # Replicate the operations in model forward pass
            src_embs = model.encoder(
                torch.tensor(all_gene_ids[i : i + batch_size], dtype=torch.long).to(
                    device
                )
            )
            val_embs = model.value_encoder(
                torch.tensor(all_values[i : i + batch_size], dtype=torch.float).to(
                    device
                )
            )
            total_embs = src_embs + val_embs
            # total_embs = model.layer(total_embs.permute(0, 2, 1)).permute(0, 2, 1)
            # Send total_embs to attention layers for attention operations
            # Retrieve the output from second to last layer
            for layer in model.transformer_encoder.layers[:num_attn_layers]:
                total_embs = layer(
                    total_embs,
                    src_key_padding_mask=src_key_padding_mask[i : i + batch_size].to(
                        device
                    ),
                )
            # Send total_embs to the last layer in flash-attn
            # https://github.com/HazyResearch/flash-attention/blob/1b18f1b7a133c20904c096b8b222a0916e1b3d37/flash_attn/flash_attention.py#L90
            try:
                qkv = F.linear(
                    total_embs,
                    model.transformer_encoder.layers[
                        num_attn_layers
                    ].self_attn.in_proj_weight,
                    bias=model.transformer_encoder.layers[
                        num_attn_layers
                    ].self_attn.in_proj_bias,
                )
            except AttributeError:
                qkv = F.linear(
                    total_embs,
                    model.transformer_encoder.layers[
                        num_attn_layers
                    ].self_attn.Wqkv.weight,
                    bias=model.transformer_encoder.layers[
                        num_attn_layers
                    ].self_attn.Wqkv.bias,
                )
            # qkv = h.reshape(total_embs.shape[0], total_embs.shape[1], nhead, 3 * d_hid//nhead)
            # qkv = qkv.permute(0, 2, 1, 3)  # [Batch, Head, SeqLen, Dims]
            # Retrieve q, k, and v from flast-attn wrapper
            qkv = rearrange(qkv, "b s (three h d) -> b s three h d", three=3, h=8)
            q = qkv[:, :, 0, :, :]
            k = qkv[:, :, 1, :, :]
            v = qkv[:, :, 2, :, :]
            # https://towardsdatascience.com/illustrated-self-attention-2d627e33b20a
            # q = [batch, gene, n_heads, n_hid]
            # k = [batch, gene, n_heads, n_hid]
            # attn_scores = [batch, n_heads, gene, gene]
            attn_scores = q.permute(0, 2, 1, 3) @ k.permute(0, 2, 3, 1)
            # apply softmax to get attention weights
            attn_scores = softmax(attn_scores, dim=-1)

            if i == 0:
                sm_attn_scores = attn_scores.sum(0).mean(0).detach().cpu().numpy()
            else:
                # take the sum
                sm_attn_scores += attn_scores.sum(0).mean(0).detach().cpu().numpy()
    sm_attn_scores = sm_attn_scores / N
    return subadata, sm_attn_scores[1:, 1:]


### prev#####


def _prepare_dataset(
    dataset,
    vocab,
    batch_size,
    epoch=5,
    mask_ratio=0.4,
    mask_value=-1,
    n_hvg=2000,
    pad_token="<pad>",
    pad_value=-2,
    per_seq_batch_sample=True,
    test_size=0.2,
    task="annotation",
):
    all_counts = (
        dataset.layers["X_binned"]
        if issparse(dataset.layers["X_binned"])
        else dataset.layers["X_binned"]
    )
    (
        train_data,
        valid_data,
        train_celltype_labels,
        valid_celltype_labels,
        train_batch_labels,
        valid_batch_labels,
    ) = train_test_split(
        all_counts,
        dataset.obs["cell_type"].astype("category").cat.codes.values,
        dataset.obs["batch_id"].astype("category").cat.codes.values,
        test_size=test_size,
        shuffle=True,
    )
    tokenized_train = tokenize_and_pad_batch(
        train_data,
        dataset.var["gene_ids"].values,
        max_len=n_hvg + 1,
        vocab=vocab,
        pad_token=pad_token,
        pad_value=pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=True,
    )
    tokenized_valid = tokenize_and_pad_batch(
        valid_data,
        dataset.var["gene_ids"].values,
        max_len=n_hvg + 1,
        vocab=vocab,
        pad_token=pad_token,
        pad_value=pad_value,
        append_cls=True,
        include_zero_gene=True,
    )
    logger.info(
        f"train set number of samples: {tokenized_train['genes'].shape[0]}, "
        f"\n\t feature length: {tokenized_train['genes'].shape[1]}"
    )
    logger.info(
        f"valid set number of samples: {tokenized_valid['genes'].shape[0]}, "
        f"\n\t feature length: {tokenized_valid['genes'].shape[1]}"
    )
    train_data_pt, valid_data_pt = prepare_data(
        tokenized_train,
        tokenized_valid,
        train_batch_labels,
        valid_batch_labels,
        task=task,
        mask_ratio=mask_ratio,
        mask_value=mask_value,
        pad_value=pad_value,
        epoch=epoch,
        train_celltype_labels=train_celltype_labels,
        valid_celltype_labels=valid_celltype_labels,
        sort_seq_batch=per_seq_batch_sample,
    )

    train_loader = prepare_dataloader(
        train_data_pt,
        batch_size=batch_size,
        shuffle=False,
        intra_domain_shuffle=True,
        drop_last=False,
    )
    valid_loader = prepare_dataloader(
        valid_data_pt,
        batch_size=batch_size,
        shuffle=False,
        intra_domain_shuffle=False,
        drop_last=False,
    )
    return train_loader, valid_loader


def setup(dataset_name, save_path, config, seed=42):
    save_dir = Path(f"{save_path}_{dataset_name}-{time.strftime('%b%d-%H-%M')}/")
    save_dir.mkdir(parents=True, exist_ok=True)
    set_seed(seed)
    add_file_handler(logger, save_dir / "run.log")
    logger.info(f"save to {save_dir}")

    return save_dir


def load_dataset(dataset, vocab):
    dataset.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in dataset.var["gene_symbols"]
    ]
    gene_ids_in_vocab = np.array(dataset.var["id_in_vocab"])
    logger.info(
        f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
        f"in vocabulary of size {len(vocab)}."
    )

    dataset = dataset[:, dataset.var["id_in_vocab"] >= 0]
    dataset.var["gene_ids"] = vocab(dataset.var["gene_symbols"].tolist())
    return dataset


def fine_tune(
    model,
    config,
    dataset,
    vocab,
    epochs,
    batch_size,
    save_folder,
    device,
    task="annotation",
):
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=config["lr"],
        eps=1e-4 if config["amp"] else 1e-8,
    )
    scheduler = torch.optim.lr_scheduler.StepLR(
        optimizer, 1, gamma=config["schedule_ratio"]
    )
    config["task"] = "annotation"
    scaler = torch.cuda.amp.GradScaler(enabled=config["amp"])

    run = wandb.init(
        config=config,
        project="scGPT",
        reinit=True,
        settings=wandb.Settings(start_method="fork"),
    )
    wandb.watch(model)

    best_val_loss = float("inf")

    best_model = None
    for epoch in range(1, epochs + 1):
        epoch_start_time = time.time()
        data_loader, valid_loader = prepare_dataset(
            dataset,
            vocab,
            batch_size,
            epoch=epoch,
            n_hvg=config["n_hvg"],
            test_size=0.2,
            mask_ratio=config["mask_ratio"],
        )
        scgpt_train(
            model,
            data_loader,
            vocab,
            criterion_gep_gepc=masked_mse_loss,
            criterion_dab=torch.nn.CrossEntropyLoss(),
            criterion_cls=nn.CrossEntropyLoss(),
            scaler=scaler,
            optimizer=optimizer,
            scheduler=scheduler,
            device=device,
            config=config,
            logger=logger,
            epoch=epoch,
        )
        val_loss = scgpt_validate(
            model,
            valid_loader,
            vocab,
            criterion_gep_gepc=masked_mse_loss,
            criterion_dab=torch.nn.CrossEntropyLoss(),
            criterion_cls=nn.CrossEntropyLoss(),
            device=device,
            config=config,
            logger=logger,
            epoch=epoch,
        )
        elapsed = time.time() - epoch_start_time
        logger.info("-" * 89)
        logger.info(
            f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
            f"valid loss/mse {val_loss:5.4f}"
        )
        logger.info("-" * 89)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model = copy.deepcopy(model)
            best_model_epoch = epoch
            logger.info(f"Best model with score {best_val_loss:5.4f}")
        scheduler.step()
    logger.info(f"Saving model to {save_folder}")

    torch.save(
        best_model.state_dict(),
        save_folder + f"model_e{best_model_epoch}.pt",
    )
    wandb.use_artifact(save_folder + f"model_e{best_model_epoch}.pt", type="model")

    wandb.finish()
    gc.collect()
    return best_model


def define_wandb_metrcis():
    wandb.define_metric("valid/mse", summary="min", step_metric="epoch")
    wandb.define_metric("valid/mre", summary="min", step_metric="epoch")
    wandb.define_metric("valid/dab", summary="min", step_metric="epoch")
    wandb.define_metric("valid/sum_mse_dab", summary="min", step_metric="epoch")
    wandb.define_metric("test/avg_bio", summary="max")


config = {
    "GEPC": False,  # Gene expression modelling for cell objective
    "MVC": True,
    "GEP": True,
    "mask_ratio": 0.4,  # Default mask ratio
    "DSBN": True,  # Default setting
    "ECS": False,
    "ecs_thres": 0.8,  # Elastic cell similarity objective, 0.0 to 1.0, 0.0 to disable
    "CLS": True,
    "USE_CLS": False,
    "USE_GENERATIVE_TRAINING": True,
    "DAR": False,
    "dab_weight": 1.0,  # DAR objective weight for batch correction
    "USE_CCE": False,
    "explicit_zero_prob": False,
    #
    "epochs": 6,  # Default number of epochs for fine-tuning
    "lr": 1e-4,  # Default learning rate for fine-tuning
    "schedule_ratio": 0.9,  # Default rate for learning rate decay
    "save_eval_interval": 5,  # Default model evaluation interval
    "log_interval": 100,  # Default log interval
    "pre_norm": False,  # Default setting
    "amp": True,  # # Default setting: Automatic Mixed Precision
    "dropout": 0.2,  # Default dropout rate during model fine-tuning
    "batch_size": 6,
    "eval_batch_size": 8,
    "log_interval": 9000,
    "save_interval": 27000,
    #
    "scheduler_interval": 100,
    "scheduler_factor": 0.99,
    "warmup_ratio_or_step": 10000.0,
    #
    "no_cce": True,
    #
    "fast_transformer": True,
    "n_layers_cls": 3,
    "fp16": True,
    "nlayers": 12,
    "embsize": 512,
    "d_hid": 512,
    "nhead": 8,  # if load model, batch_size, layer_size, nlayers, nhead will be ignored
    "layer_size": 128,
    #
    "max_seq_len": 1200,
    "n_hvg": 2000,
    "n_bins": 51,  # Default number of bins for value binning in data pre-processing
    "mask_value": -1,
    "pad_value": -2,
    "pad_token": "<pad>",
    "input_style": "binned",
    #
    "valid_size_or_ratio": 0.003,
    "dist_backend": "nccl",
    "grad_accu_steps": 1,
    "input_emb_style": "continuous",
    "training_tasks": "both",
    "trunc_by_sample": True,
    "rank": 0,
    #
    "world_size": 16,
    "distributed": True,
    "local_rank": 0,
    "gpu": 0,
    "task": "annotation",
    "use_batch_labels": True,
    "use_mod": False,
}
