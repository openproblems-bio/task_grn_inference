#!/usr/bin/env python3
"""
Filter TF binding datasets to keep only specific cell types.
Supports multiple cell type categories including PBMC, K562, HEK293, HCT116, etc.
"""

import pandas as pd
import json
import re
import os
import argparse
from pathlib import Path

def load_cell_type_categories():
    """Load cell type categories for different experimental systems"""
    
    cell_type_cats = {
        # PBMC-related cell types (immune cells)
        "pbmc": {
            "tfmdb": ["B cell", "Dendritic cell", "Lymphoid cell", "Macrophage cell", "Macrophages cell", "Macrophagocyte cell", "Mononuclear cell", "Natural killer cell", "Peripheral blood cell", "T cell"],
            "hpa": ["B", "Dendritic cells", "Erythroid cells", "Lymphoid tissue", "Macrophages", "Monocytes", "NK", "T"],
            "chipatlas": ["B-Cell", "B-cell (CD19+)", "Bone Marrow Cells", "Bone marrow mononuclear cells", "Bone marrow nuclear cells", "CD20+ B cells", "CD4+ T cells", "CD8+ T cells", 
                        "Dendritic Cells", "Macrophages", "Mast Cells", "Memory B cells", "Memory T cells", "Monocytes", "Monocytes-CD14+", "Naive B cells", "Natural Killer T-Cells", 
                        "Natural killer cells", "T cells", "Th1 Cells", "Th17 Cells", "Th2 Cells", "Treg", "Blood", "GM12878"],  # Added Blood and GM12878
            "remap2022": ["B-cell", "CD4", "CD4-pos", "CD8", "T-cell", "Th1", "Th17", "macrophage", "monocyte", "peripheral-blood-mononuclear-cell", "primary-B-cell", "primary-monocyte"],
            "unibind": ["B cells", "CD4 CD25 CD45RA T cells", "CD4 CD25 T cells", "CD4 T cells", "T cells", "Th17 cells", "Th2 cells", "macrophage", "macrophage hypo il", "macrophage il4", 
                    "macrophages", "macrophages ifng", "macrophages il4", "macrophages tpp", "monocyte derived macrophages", "monocyte macro", "monocytes blood", "monocytes of peripheral blood", 
                    "peripheral blood mononuclear cells", "GM12878"],  # Added GM12878
        },
        

        # K562 cell line (myelogenous leukemia)
        "k562": {
            "chipatlas": ["K-562"],
            "remap2022": ["K-562 cell"],
            "unibind": ["K562", "K562  myelogenous leukemia"],
        },
        
        # HEK293 cell lines (embryonic kidney)
        "hek293": {
            "chipatlas": ["293", "HEK293"],  # Found "293,Kidney" in data
            "remap2022": ["HEK-293 cell", "HEK-293T cell"],
            "unibind": ["HEK293", "HEK293  embryonic kidney"],
        },
        
        # HCT116 cell line (colorectal carcinoma)
        "hct116": {
            "chipatlas": ["HCT116", "HCT-116", "DLD-1", "HT-29"],
            "remap2022": ["HCT-116 cell", "DLD-1 cell", "CACO-2 cell"],
            "unibind": ["HCT116"],
        }
    }
    
    return cell_type_cats

def filter_dataset(input_file, output_file, cell_terms, database_name, category_name):
    """Filter a TF binding dataset for specific cell types"""
    print(f"Filtering {database_name} dataset for {category_name}...")
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Cell terms: {len(cell_terms)} terms")
    print(f"Search terms: {cell_terms}")
    
    # Read the dataset
    df = pd.read_csv(input_file, sep='\t', header=None, 
                     names=['Chromosome', 'Start', 'End', 'TF', 'CellType'])
    
    print(f"Original dataset: {len(df)} entries")
    
    def has_target_cell_type(cell_string):
        """Check if entry contains ANY target cell types"""
        if pd.isna(cell_string):
            return False
        
        # Split by comma and strip whitespace
        individual_cells = [cell.strip() for cell in cell_string.split(',')]
        
        # Check if ANY individual cell type matches our target terms
        for cell in individual_cells:
            cell_lower = cell.lower()
            for term in cell_terms:
                term_lower = term.lower()
                # Use exact match or contains match depending on specificity
                if term_lower == cell_lower or term_lower in cell_lower:
                    return True
        
        return False
    
    # Apply filtering
    filtered_df = df[df['CellType'].apply(has_target_cell_type)]
    
    print(f"After {category_name} filtering: {len(filtered_df)} entries ({len(filtered_df)/len(df)*100:.1f}%)")
    
    # Count unique TFs
    original_tfs = df['TF'].nunique()
    filtered_tfs = filtered_df['TF'].nunique()
    print(f"Original TFs: {original_tfs}")
    print(f"Filtered TFs: {filtered_tfs} ({filtered_tfs/original_tfs*100:.1f}%)")
    
    # Save filtered dataset
    if len(filtered_df) > 0:
        filtered_df.to_csv(output_file, sep='\t', header=False, index=False)
        print(f"Saved {len(filtered_df)} entries to {output_file}")
    else:
        print(f"No entries found for {category_name} in {database_name}")
        return
    
    # Show some examples of kept cell types
    unique_cells = filtered_df['CellType'].unique()[:10]
    print(f"Example kept cell types: {list(unique_cells)}")
    
    # Show TF distribution
    tf_counts = filtered_df['TF'].value_counts().head(10)
    print(f"Top TFs in {category_name} data: {dict(tf_counts)}")
    print()

def main():
    """Main function to filter TF binding datasets for specific cell types"""
    
    parser = argparse.ArgumentParser(description='Filter TF binding datasets by cell type')
    parser.add_argument('--cell_type', type=str, default='pbmc',
                       choices=['pbmc', 'k562', 'hek293', 'hct116', 'mcf7', 'hepg2', 'a549', 'cancer_lines'],
                       help='Cell type category to filter for')
    parser.add_argument('--all', action='store_true',
                       help='Generate filtered datasets for all cell type categories')
    
    args = parser.parse_args()
    
    # Load cell type categories
    cell_type_cats = load_cell_type_categories()
    
    # Define base paths
    base_input_dir = 'resources/datasets_raw/tfb'
    base_output_dir = 'resources/datasets_raw/tfb'
    
    datasets = ['remap2022', 'chipatlas', 'unibind']
    
    # Determine which categories to process
    if args.all:
        categories_to_process = list(cell_type_cats.keys())
    else:
        categories_to_process = [args.cell_type]
    
    print(f"Filtering TF binding datasets for: {categories_to_process}")
    print("=" * 60)
    
    # Process each category
    for category in categories_to_process:
        print(f"\n=== PROCESSING {category.upper()} ===")
        
        if category not in cell_type_cats:
            print(f"Unknown category: {category}")
            continue
            
        category_terms = cell_type_cats[category]
        
        # Filter each dataset
        for db_name in datasets:
            input_file = f"{base_input_dir}/{db_name}/{db_name}.bed"
            output_file = f"{base_output_dir}/{db_name}/{db_name}_{category}.bed"
            
            if os.path.exists(input_file):
                # Create output directory if needed
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                
                # Get cell terms for this database
                if db_name in category_terms:
                    terms = category_terms[db_name]
                    
                    filter_dataset(
                        input_file=input_file,
                        output_file=output_file,
                        cell_terms=terms,
                        database_name=db_name,
                        category_name=category
                    )
                else:
                    print(f"No {category} terms defined for {db_name}, skipping...")
            else:
                print(f"Warning: {input_file} not found, skipping {db_name}")
    
    print("\nFiltering complete!")
    print(f"\nFiltered datasets created in: {base_output_dir}")
    
    # Show summary of created files
    for category in categories_to_process:
        if category in cell_type_cats:
            print(f"\n{category.upper()} datasets:")
            for db_name in datasets:
                output_file = f"{base_output_dir}/{db_name}/{db_name}_{category}.bed"
                if os.path.exists(output_file):
                    # Count lines
                    with open(output_file, 'r') as f:
                        line_count = sum(1 for _ in f)
                    print(f"  - {output_file} ({line_count:,} entries)")

if __name__ == '__main__':
    main()