#!/usr/bin/env python3
"""
Download and save pathway gene sets from Enrichr for GRN evaluation.

This script downloads pathway databases and saves them as GMT files
in the resources directory for later use in annotation metrics.

Usage:
    python acquisition.py --output_dir resources/grn_benchmark/prior/pathways/
"""

import os
import sys
import argparse
import requests
from pathlib import Path
from typing import Dict, List
from tqdm import tqdm


def get_enrichr_library(library_name: str) -> Dict[str, List[str]]:
    """
    Download gene sets from Enrichr.
    
    Parameters
    ----------
    library_name : str
        Name of the Enrichr library
    
    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping pathway names to gene lists
    """
    url = f'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}'
    response = requests.get(url)
    
    if response.status_code != 200:
        raise ValueError(f"Failed to download {library_name}: {response.status_code}")
    
    gene_sets = {}
    for line in response.text.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        term = parts[0]
        genes = [g for g in parts[2:] if g]  # Filter empty strings
        if genes:  # Only add if there are genes
            gene_sets[term] = genes
    
    return gene_sets


def save_as_gmt(gene_sets: Dict[str, List[str]], output_file: str):
    """
    Save gene sets in GMT format.
    
    Parameters
    ----------
    gene_sets : Dict[str, List[str]]
        Gene sets to save
    output_file : str
        Output GMT file path
    """
    with open(output_file, 'w') as f:
        for term, genes in gene_sets.items():
            # GMT format: name\tdescription\tgene1\tgene2\t...
            line = f"{term}\t{term}\t" + "\t".join(genes) + "\n"
            f.write(line)


def get_library_stats(gene_sets: Dict[str, List[str]], 
                     min_size: int = 10, 
                     max_size: int = 500) -> Dict:
    """
    Calculate statistics for a gene set library.
    
    Parameters
    ----------
    gene_sets : Dict[str, List[str]]
        Gene sets
    min_size : int
        Minimum pathway size
    max_size : int
        Maximum pathway size
    
    Returns
    -------
    Dict
        Statistics including total, filtered counts, size distribution
    """
    sizes = [len(genes) for genes in gene_sets.values()]
    filtered = [s for s in sizes if min_size <= s <= max_size]
    
    return {
        'total_pathways': len(gene_sets),
        'filtered_pathways': len(filtered),
        'min_size': min(sizes) if sizes else 0,
        'max_size': max(sizes) if sizes else 0,
        'median_size': sorted(sizes)[len(sizes)//2] if sizes else 0,
        'mean_size': sum(sizes) / len(sizes) if sizes else 0
    }


def main():
    parser = argparse.ArgumentParser(
        description='Download pathway gene sets from Enrichr'
    )
    parser.add_argument(
        '--output_dir', 
        default='resources/grn_benchmark/prior/pathways/',
        help='Output directory for GMT files'
    )
    parser.add_argument(
        '--min_size', 
        type=int, 
        default=10, 
        help='Minimum pathway size for filtering statistics'
    )
    parser.add_argument(
        '--max_size', 
        type=int, 
        default=500, 
        help='Maximum pathway size for filtering statistics'
    )
    
    args = parser.parse_args()
    
    # Pathway databases to download
    libraries = {
        'MSigDB_Hallmark_2020': 'hallmark_2020',
        'GO_Biological_Process_2023': 'go_bp_2023',
        'KEGG_2021_Human': 'kegg_2021',
        'WikiPathways_2019_Human': 'wikipathways_2019',
        'BioPlanet_2019': 'bioplanet_2019',
        'Reactome_2022': 'reactome_2022',
    }
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("PATHWAY GENE SET ACQUISITION")
    print("=" * 80)
    print(f"\nOutput directory: {output_dir}")
    print(f"Downloading {len(libraries)} pathway databases...\n")
    
    # Download and save each library
    summary = []
    
    for library_name, short_name in tqdm(libraries.items(), desc="Downloading libraries"):
        print(f"\n{'='*80}")
        print(f"Processing: {library_name}")
        print(f"{'='*80}")
        
        try:
            # Download
            print(f"  [1/3] Downloading from Enrichr...")
            gene_sets = get_enrichr_library(library_name)
            print(f"    Downloaded {len(gene_sets)} gene sets")
            
            # Calculate statistics
            print(f"  [2/3] Calculating statistics...")
            stats = get_library_stats(gene_sets, args.min_size, args.max_size)
            print(f"    Total pathways: {stats['total_pathways']}")
            print(f"    After filtering ({args.min_size}-{args.max_size}): {stats['filtered_pathways']}")
            print(f"    Size range: {stats['min_size']}-{stats['max_size']} (median: {stats['median_size']:.0f})")
            
            # Save as GMT
            output_file = output_dir / f"{short_name}.gmt"
            print(f"  [3/3] Saving to {output_file}...")
            save_as_gmt(gene_sets, str(output_file))
            print(f"    Saved successfully")
            
            # Add to summary
            summary.append({
                'library': library_name,
                'short_name': short_name,
                'file': str(output_file),
                **stats
            })
            
        except Exception as e:
            print(f"    ERROR: {e}")
            summary.append({
                'library': library_name,
                'short_name': short_name,
                'file': 'FAILED',
                'error': str(e)
            })
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    for item in summary:
        if 'error' in item:
            print(f"\n{item['library']} ({item['short_name']}): FAILED")
            print(f"  Error: {item['error']}")
        else:
            print(f"\n{item['library']} ({item['short_name']}):")
            print(f"  File: {item['file']}")
            print(f"  Total pathways: {item['total_pathways']}")
            print(f"  Filtered ({args.min_size}-{args.max_size}): {item['filtered_pathways']}")
            print(f"  Size range: {item['min_size']}-{item['max_size']} (median: {item['median_size']:.0f})")
    
    print("\n" + "=" * 80)
    print("ACQUISITION COMPLETE")
    print("=" * 80)
    print(f"\nPathway files saved in: {output_dir}")
    print("\nFiles created:")
    for item in summary:
        if 'error' not in item:
            print(f"  - {item['short_name']}.gmt ({item['total_pathways']} pathways)")
    
    print("\n" + "=" * 80)


if __name__ == '__main__':
    main()
