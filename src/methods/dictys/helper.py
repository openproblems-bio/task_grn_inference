import os
os.environ["MKL_SERVICE_FORCE_INTEL"] = "1"
os.environ["MKL_THREADING_LAYER"] = "GNU"
import shutil
from typing import Optional, List

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import warnings
import subprocess
warnings.filterwarnings("ignore")


OVERRIDE_DOWNLOAD = True


def run_cmd(cmd: List[str], cwd: Optional[str] = None) -> None:
    kwargs = {}
    if cwd is not None:
        kwargs['cwd'] = cwd
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        **kwargs
    ) as proc:
        for line in proc.stdout:
            print(line, end="")
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"Command {cmd} failed with exit code {rc}")


def define_vars(par):
    os.makedirs(par['temp_dir'], exist_ok=True)
    print(f"Temporary directory: {par['temp_dir']}", flush=True)
    par['data_dir'] = f"{par['temp_dir']}/data/"
    os.makedirs(par['data_dir'], exist_ok=True)
    par['exp_path'] = f"{par['data_dir']}/expression.tsv.gz"
    par['barcodes'] = f"{par['data_dir']}/barcodes_all.txt"
    par['pks_path'] = f"{par['data_dir']}/peaks.bed"
    par['frags_path'] = f"{par['data_dir']}/fragments.tsv.gz"

    par['bam_name'] = f"{par['data_dir']}/reads_all.bam"
    par['bai_name'] = f"{par['data_dir']}/reads_all.bai"

    par['bams_dir'] = f"{par['data_dir']}/bams/"

    par['gene_bed'] = f"{par['data_dir']}/gene.bed"
    par['make_dir'] = f"{par['temp_dir']}/makefiles"
    

def extract_exp(par):
    print('Extracting expression matrix and barcodes', flush=True)
    import anndata as ad
    rna = ad.read_h5ad(par['rna'])
    if 'counts' not in rna.layers:
        X = rna.X
    else:
        X = rna.layers['counts']
    if not isinstance(X, np.ndarray):
        X = X.todense()
    barcodes = rna.obs.index
    # barcodes = [f"CB:Z:{bc.replace('-', '')} for bc in barcodes]
    rna_X = pd.DataFrame(X.T, columns=barcodes, index=rna.var.index)
    rna_X.to_csv(par['exp_path'], sep="\t", compression="gzip")
    all_ids = barcodes
    with open(par['barcodes'], "w") as f:
        for i in all_ids:
            f.write(f"{i}\n")

def extract_atac(par):
    print('Extracting peaks and fragments', flush=True)
    import anndata as ad
    atac = ad.read_h5ad(par['atac'])
    all_atac_peak = np.unique([n.replace(':', '-') for n in atac.var.index])
    all_atac_peak = pd.DataFrame([n.split('-') for n in all_atac_peak])
    all_atac_peak.columns = ['chr', 'srt', 'end']
    all_atac_peak['srt'] = all_atac_peak['srt'].astype(int)
    all_atac_peak['end'] = all_atac_peak['end'].astype(int)
    all_atac_peak = all_atac_peak[(all_atac_peak.end - all_atac_peak.srt) >= 100]
    all_atac_peak = all_atac_peak.sort_values(by=['chr', 'srt', 'end'])
    all_atac_peak.to_csv(par['pks_path'], sep='\t', header=False, index=False)
    # get fragments
    coo_matrix = atac.X.tocoo()
    rows = coo_matrix.row
    cols = coo_matrix.col
    counts = coo_matrix.data
    row_names = atac.obs.index[rows]
    coord_names = atac.var.index[cols]

    # Build DataFrame
    d = {
        'chromosome': [s.split(':')[0] for s in coord_names],
        'start': [int(s.split(':')[1].split('-')[0]) for s in coord_names],
        'end': [int(s.split(':')[1].split('-')[1]) for s in coord_names],
        'barcode': [barcode.strip('\n') for barcode in row_names],
        'count': np.squeeze(np.asarray(counts.astype(np.uint16)))
    }
    df = pd.DataFrame(d, copy=False)

    # Write to TSV without header/index
    frags_path = par['frags_path']
    temp_path = frags_path.rstrip('.gz')
    df.to_csv(temp_path, sep='\t', index=False, header=False)

    # Compress TSV
    print(f'Sort and compress tsv file {frags_path}')
    os.system(f"sort -k1,1 -k2,2n {temp_path} | bgzip -c > {frags_path}")


def create_bam(par):
    print('Creating BAM file from fragments', flush=True)
    cmd = f"python {par['frag_to_bam']} --fnames {par['frags_path']} --barcodes {par['barcodes']}"
    view_sort_cmd = f"samtools view -b | samtools sort -o {par['bam_name']}"
    full_cmd = f"{cmd} | {view_sort_cmd}"
    subprocess.run(full_cmd, shell=True, check=True)
    subprocess.run(["samtools", "index", par['bam_name'], par['bai_name']], check=True)
    

def bam_to_bams(par):
    """
    Split a BAM file into per-cell BAMs using the split_bam.sh script.

    Args:
        par: dict with keys:
            - 'bam_name': path to the input BAM file
            - 'bams_dir': path to output folder for per-cell BAMs
            - 'exp_path': path to reference expression file
    """

    print('Delete temp BAM directories', flush=True)
    folders = [
        par['bams_dir'],
        os.path.join(par['bams_dir'], '..', 'bams_text'),
        os.path.join(par['bams_dir'], '..', 'bams_header')
    ]
    for folder in folders:
        if os.path.exists(folder):
            shutil.rmtree(folder)

    print('Splitting BAM into per-cell BAMs', flush=True)
    run_cmd([
        "bash", "dictys_helper", "split_bam.sh", par['bam_name'], par['bams_dir'],
        "--section", "CB:Z:", "--ref_expression", par['exp_path']
    ])


def extrac_clusters(par):
    print('Extracting clusters', flush=True)
    subsets = f"{par['data_dir']}/subsets.txt"
    with open(subsets, "w") as f:
        f.write("cluster_1\n")
    os.makedirs(f"{par['data_dir']}/subsets/", exist_ok=True)
    subsets_sub = f"{par['data_dir']}/subsets/cluster_1/"
    os.makedirs(subsets_sub, exist_ok=True)
    barcode_rna = f"{subsets_sub}/names_rna.txt"
    barcode_atac = f"{subsets_sub}/names_atac.txt"

    cp = f"cp {par['barcodes']} {barcode_rna}"
    subprocess.run(cp, shell=True, check=True)
    cp = f"cp {par['barcodes']} {barcode_atac}"
    subprocess.run(cp, shell=True, check=True)
    print('Extracting clusters successful', flush=True)


def download_file(url, dest):
    import requests
    print(f"Downloading {url} to {dest}", flush=True)
    with requests.get(url, stream=True) as r:
        r.raise_for_status()  # will raise an error if the download fails
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


def get_priors(par):
    import gzip
    import shutil
    # - get the genome
    print('Getting genome ...', flush=True)
    if OVERRIDE_DOWNLOAD or (not os.path.exists(f"{par['data_dir']}/genome/genome.fa")):
        os.makedirs(f"{par['data_dir']}/genome/", exist_ok=True)
        try:
            run_cmd([
                "aws", "s3", "cp", "s3://openproblems-data/resources/grn/supp_data/genome/genome.fa",
                f"{par['data_dir']}/genome/", "--no-sign-request"
            ])
        except:
            try:
                run_cmd([
                    "cp", "resources/supp_data/genome/genome.fa", f"{par['data_dir']}/genome/"
                ])
            except:
                raise ValueError("Could not get the reference genome")
    
    # - get gene annotation       
    print('Getting gene annotation ...', flush=True)
    data_dir = Path(par['data_dir'])
    gtf_gz = data_dir / "gene.gtf.gz"
    gtf =  data_dir / "gene.gtf"
    url = "http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz"
    if OVERRIDE_DOWNLOAD or (not os.path.exists(gtf_gz)):
        download_file(url, gtf_gz)

    with gzip.open(gtf_gz, "rb") as f_in, open(gtf, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gtf_gz.unlink()

    print('Making bed files for gene annotation ...', flush=True)
    run_cmd([
        "bash", "dictys_helper", "gene_gtf.sh", gtf, par['gene_bed']
    ])

    print('Downloading motif file...', flush=True)

    url='https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif'
    motif_file = data_dir / 'motifs.motif'
    if OVERRIDE_DOWNLOAD or (not os.path.exists(motif_file)):
        download_file(url, motif_file)


def configure(par):
    import json
    device='cuda:0' #cuda:0 , cpu
    os.makedirs(par['make_dir'], exist_ok=True)
    run_cmd([
        "bash", "dictys_helper", "makefile_template.sh", "common.mk", "config.mk", "env_none.mk", "static.mk"
    ], cwd=par['make_dir'])
    
    json_arg = json.dumps({
        "DEVICE": device,
        "GENOME_MACS2": "hs",
        "JOINT": "1"
    })

    run_cmd([
        "bash", "dictys_helper", "makefile_update.py", "config.mk", json_arg
    ], cwd=par['make_dir'])

    run_cmd([
        "bash", "dictys_helper", "makefile_check.py", "--dir_makefiles", par['make_dir'],
        "--dir_data", par['data_dir']
    ])


def infer_grn(par):
    print('Inferring GRNs', flush=True)
    run_cmd([
        "bash", "dictys_helper", "network_inference.sh", "-j", str(par['num_workers']), "-J", "1", "static"
    ], cwd=par['temp_dir'])


def export_net(par):
    from util import process_links
    from dictys.net import network
    d0 = network.from_file(f"{par['temp_dir']}/output/static.h5") # check how to control this
    d0.export(f"{par['temp_dir']}/output/net",sparsities=[None])  # no way to control the name it's exporting
    net = pd.read_csv(f"{par['temp_dir']}/output/net/Full/cluster_1.tsv.gz", sep='\t') # this is customized

    prediction = net.reset_index().melt(id_vars=net.index.name or 'index', 
                                 var_name='target', 
                                 value_name='weight')
    net.rename(columns={net.index.name or 'index': 'source'}, inplace=True)

    net = process_links(net)

    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": "dictys", "dataset_id": ad.read_h5ad(par['rna'], backed='r').uns['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])

def main(par):

    define_vars(par)
    extract_exp(par)
    extract_atac(par)
    create_bam(par)
    bam_to_bams(par)
    extrac_clusters(par)
    get_priors(par)
    configure(par)
    infer_grn(par)
    export_net(par)
