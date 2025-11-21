"""Utilities to download gene set libraries from Enrichr and save as GMT.

This is intentionally lightweight and dependency-free apart from `requests`.
"""
from typing import Dict, List
import requests


def get_enrichr_library(library_name: str) -> Dict[str, List[str]]:
    """Download a gene set library from Enrichr.

    Parameters
    ----------
    library_name: str
        Name of the Enrichr library (e.g. 'GO_Biological_Process_2023')

    Returns
    -------
    dict
        Mapping from term name -> list of gene symbols
    """
    url = f'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}'
    resp = requests.get(url)
    resp.raise_for_status()

    lines = resp.text.strip().split('\n')
    gene_sets = {}
    for line in lines:
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        term = parts[0]
        genes = [g for g in parts[2:] if g]
        gene_sets[term] = genes

    return gene_sets


def save_library_gmt(gene_sets: Dict[str, List[str]], out_path: str) -> None:
    """Save a gene sets dict to GMT format.

    Parameters
    ----------
    gene_sets: dict
        Mapping term -> list of genes
    out_path: str
        Path to write GMT file
    """
    with open(out_path, 'w') as f:
        for term, genes in gene_sets.items():
            # GMT: term <tab> description <tab> gene1 <tab> gene2 ...
            desc = term
            fields = [term, desc] + genes
            f.write('\t'.join(fields) + '\n')


if __name__ == '__main__':
    # simple CLI for quick manual use
    import sys
    if len(sys.argv) < 3:
        print('Usage: python enrichr_utils.py <library_name> <out_gmt>')
        sys.exit(1)
    lib = sys.argv[1]
    out = sys.argv[2]
    gs = get_enrichr_library(lib)
    save_library_gmt(gs, out)
    print(f'Wrote {len(gs)} gene sets to {out}')
