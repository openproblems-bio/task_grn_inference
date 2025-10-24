#!/usr/bin/env python3
"""
Create two versions of merged PBMC ground truth networks:
1. Union: All edges from any dataset (what we had before)
2. Intersection: Only edges common to ALL datasets
"""

import pandas as pd
import numpy as np

def create_union_and_intersection_gt():
    """Create union and intersection GT networks from the three PBMC datasets."""
    
    base_dir = "resources/grn_benchmark/ground_truth/"
    
    # Load the three GT datasets
    print("=== Loading PBMC GT datasets ===")
    datasets = {}
    
    for name, file in [('UniBind', 'PBMC_unibind_pbmc.csv'), 
                       ('ChIP-Atlas', 'PBMC_chipatlas_pbmc.csv'), 
                       ('ReMap2022', 'PBMC_remap2022_pbmc.csv')]:
        try:
            df = pd.read_csv(f"{base_dir}/{file}")
            datasets[name] = df
            print(f"{name:12}: {len(df):6,} edges, {df['source'].nunique():2d} TFs, {df['target'].nunique():5,} targets")
        except FileNotFoundError:
            print(f"Warning: {file} not found")
            continue
    
    if len(datasets) != 3:
        print("Need all 3 datasets for intersection analysis!")
        return
    
    # Convert to edge sets for set operations
    print(f"\n=== Converting to edge sets ===")
    edge_sets = {}
    
    for name, df in datasets.items():
        edges = set(zip(df['source'], df['target']))
        edge_sets[name] = edges
        print(f"{name:12}: {len(edges):6,} unique edges")
    
    # Create UNION (all edges from any dataset)
    print(f"\n=== CREATING UNION GT ===")
    union_edges = set()
    for edges in edge_sets.values():
        union_edges.update(edges)
    
    print(f"Union edges: {len(union_edges):,}")
    
    # Convert union back to DataFrame
    union_df = pd.DataFrame(list(union_edges), columns=['source', 'target'])
    union_df = union_df.sort_values(['source', 'target']).reset_index(drop=True)
    
    # Union statistics
    n_tfs_union = union_df['source'].nunique()
    n_targets_union = union_df['target'].nunique()
    n_edges_union = len(union_df)
    
    print(f"Union network:")
    print(f"  TFs:     {n_tfs_union:2d}")
    print(f"  Targets: {n_targets_union:5,}")
    print(f"  Edges:   {n_edges_union:6,}")
    print(f"  Density: {n_edges_union/(n_tfs_union*n_targets_union)*100:.2f}%")
    
    # Create INTERSECTION (only edges common to ALL datasets)
    print(f"\n=== CREATING INTERSECTION GT ===")
    intersection_edges = set.intersection(*edge_sets.values())
    
    print(f"Intersection edges: {len(intersection_edges):,}")
    
    # Convert intersection back to DataFrame
    intersection_df = pd.DataFrame(list(intersection_edges), columns=['source', 'target'])
    intersection_df = intersection_df.sort_values(['source', 'target']).reset_index(drop=True)
    
    # Intersection statistics
    n_tfs_intersect = intersection_df['source'].nunique()
    n_targets_intersect = intersection_df['target'].nunique()
    n_edges_intersect = len(intersection_df)
    
    print(f"Intersection network:")
    print(f"  TFs:     {n_tfs_intersect:2d}")
    print(f"  Targets: {n_targets_intersect:5,}")
    print(f"  Edges:   {n_edges_intersect:6,}")
    if n_tfs_intersect > 0 and n_targets_intersect > 0:
        print(f"  Density: {n_edges_intersect/(n_tfs_intersect*n_targets_intersect)*100:.2f}%")
    
    # Analyze which TFs are in intersection
    print(f"\n=== INTERSECTION TF ANALYSIS ===")
    intersect_tfs = sorted(intersection_df['source'].unique())
    print(f"TFs in intersection ({len(intersect_tfs)}): {intersect_tfs}")
    
    # Compare TF coverage
    union_tfs = set(union_df['source'].unique())
    print(f"\nTF coverage comparison:")
    print(f"  Union TFs:        {len(union_tfs):2d}")
    print(f"  Intersection TFs: {len(intersect_tfs):2d}")
    print(f"  TFs lost:         {len(union_tfs) - len(intersect_tfs):2d}")
    
    # Show degree distributions
    print(f"\n=== DEGREE DISTRIBUTIONS ===")
    
    # Union degrees
    union_tf_degrees = union_df.groupby('source').size().sort_values(ascending=False)
    union_target_degrees = union_df.groupby('target').size().sort_values(ascending=False)
    
    print(f"Union network:")
    print(f"  TF degrees    - Min: {union_tf_degrees.min():4,}, Max: {union_tf_degrees.max():5,}, Mean: {union_tf_degrees.mean():6.1f}")
    print(f"  Target degrees - Min: {union_target_degrees.min():2d}, Max: {union_target_degrees.max():2d}, Mean: {union_target_degrees.mean():4.1f}")
    
    # Intersection degrees (if not empty)
    if not intersection_df.empty:
        intersect_tf_degrees = intersection_df.groupby('source').size().sort_values(ascending=False)
        intersect_target_degrees = intersection_df.groupby('target').size().sort_values(ascending=False)
        
        print(f"Intersection network:")
        print(f"  TF degrees    - Min: {intersect_tf_degrees.min():4,}, Max: {intersect_tf_degrees.max():5,}, Mean: {intersect_tf_degrees.mean():6.1f}")
        print(f"  Target degrees - Min: {intersect_target_degrees.min():2d}, Max: {intersect_target_degrees.max():2d}, Mean: {intersect_target_degrees.mean():4.1f}")
        
        print(f"\nTop 10 TFs in intersection by degree:")
        for i, (tf, degree) in enumerate(intersect_tf_degrees.head(10).items(), 1):
            print(f"  {i:2d}. {tf:<8}: {degree:4,} targets")
    
    # Analyze what's lost in intersection
    print(f"\n=== WHAT'S LOST IN INTERSECTION ===")
    
    union_only_edges = union_edges - intersection_edges
    print(f"Edges lost: {len(union_only_edges):,} ({len(union_only_edges)/len(union_edges)*100:.1f}% of union)")
    
    # Which TFs lose the most edges?
    union_only_df = pd.DataFrame(list(union_only_edges), columns=['source', 'target'])
    if not union_only_df.empty:
        lost_by_tf = union_only_df.groupby('source').size().sort_values(ascending=False)
        print(f"\nTFs losing most edges in intersection:")
        for i, (tf, lost_count) in enumerate(lost_by_tf.head(10).items(), 1):
            union_count = union_tf_degrees.get(tf, 0)
            intersect_count = intersect_tf_degrees.get(tf, 0) if not intersection_df.empty else 0
            pct_lost = lost_count / union_count * 100 if union_count > 0 else 0
            print(f"  {i:2d}. {tf:<8}: {lost_count:4,} edges lost ({pct_lost:4.1f}% of {union_count:,})")
    
    # Save both networks
    print(f"\n=== SAVING NETWORKS ===")
    
    # Union network
    union_file = f"{base_dir}/PBMC_union_all3.csv"
    union_df.to_csv(union_file, index=False)
    print(f"Union GT saved to:        {union_file}")
    print(f"  Format: CSV with 'source', 'target' columns")
    print(f"  Size: {len(union_df):,} edges")
    
    # Intersection network
    intersect_file = f"{base_dir}/PBMC_intersection_all3.csv"
    intersection_df.to_csv(intersect_file, index=False)
    print(f"Intersection GT saved to: {intersect_file}")
    print(f"  Format: CSV with 'source', 'target' columns")
    print(f"  Size: {len(intersection_df):,} edges")
    
    # Summary comparison
    print(f"\n=== FINAL COMPARISON ===")
    print(f"{'Network':<12} {'TFs':<4} {'Targets':<8} {'Edges':<8} {'Density'}")
    print(f"{'-'*50}")
    
    for name, df in datasets.items():
        density = len(df) / (df['source'].nunique() * df['target'].nunique()) * 100
        print(f"{name:<12} {df['source'].nunique():2d}   {df['target'].nunique():7,}  {len(df):7,}  {density:6.2f}%")
    
    union_density = n_edges_union / (n_tfs_union * n_targets_union) * 100
    print(f"{'Union':<12} {n_tfs_union:2d}   {n_targets_union:7,}  {n_edges_union:7,}  {union_density:6.2f}%")
    
    if n_tfs_intersect > 0 and n_targets_intersect > 0:
        intersect_density = n_edges_intersect / (n_tfs_intersect * n_targets_intersect) * 100
        print(f"{'Intersection':<12} {n_tfs_intersect:2d}   {n_targets_intersect:7,}  {n_edges_intersect:7,}  {intersect_density:6.2f}%")
    
    return union_df, intersection_df

if __name__ == '__main__':
    union_gt, intersection_gt = create_union_and_intersection_gt()