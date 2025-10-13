import pandas as pd

for cell_type in ['K562']: #'K562'
    print(f'Processing cell type: {cell_type}', flush=True)

    gt_1 = pd.read_csv(f"resources/grn_benchmark/ground_truth/{cell_type}_chipatlas.csv")
    gt_2 = pd.read_csv(f"resources/grn_benchmark/ground_truth/{cell_type}_remap.csv")

    # Ensure both DataFrames contain the required columns
    for col in ['source', 'target']:
        if col not in gt_1.columns or col not in gt_2.columns:
            raise ValueError(f"Missing column '{col}' in one of the dataframes")

    # Convert to strings to avoid NaN / mixed-type issues
    gt_1['source'] = gt_1['source'].astype(str)
    gt_1['target'] = gt_1['target'].astype(str)
    gt_2['source'] = gt_2['source'].astype(str)
    gt_2['target'] = gt_2['target'].astype(str)

    # Compute overlaps
    source_overlap = set(gt_1['source']).intersection(gt_2['source'])
    target_overlap = set(gt_1['target']).intersection(gt_2['target'])

    # Source-target pair overlap
    gt_1_pairs = set(zip(gt_1['source'], gt_1['target']))
    gt_2_pairs = set(zip(gt_2['source'], gt_2['target']))
    pair_overlap = gt_1_pairs.intersection(gt_2_pairs)

    # Display results
    print(f"Overlapping sources: {len(source_overlap)} / {len(set(gt_1['source']).union(gt_2['source']))}")
    print(f"Overlapping targets: {len(target_overlap)} / {len(set(gt_1['target']).union(gt_2['target']))}")
    print(f"Overlapping source-target pairs: {len(pair_overlap)} / {len(set(gt_1_pairs).union(gt_2_pairs))}")

    # Save overlapping edges
    mutual_df = pd.DataFrame(list(pair_overlap), columns=['source', 'target'])
    mutual_df.to_csv(f"resources/grn_benchmark/ground_truth/{cell_type}_mutual.csv", index=False)

    print(f"Saved {len(mutual_df)} overlapping edges to 'resources/grn_benchmark/ground_truth/{cell_type}_mutual.csv'")