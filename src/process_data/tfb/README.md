# TFB Dataset Processing

This directory contains scripts to process TF binding (TFB) datasets and create Ground Truth Gene Regulatory Networks (GRNs) for multiple cell types.

## Files

- `process_tfb_datasets.py` - Main script to process TFB datasets for multiple cell types
- `build_grn.py` - Core functions for building GRNs from TFB data
- `README.md` - This file

## Usage

### List Available Cell Types

```bash
python process_tfb_datasets.py --list-available
```

### Process Specific Cell Types

Process a single cell type:
```bash
python process_tfb_datasets.py --cell-types k562
```

Process multiple cell types:
```bash
python process_tfb_datasets.py --cell-types pbmc k562 hek293
```

### Process All Available Cell Types

```bash
python process_tfb_datasets.py
```

## Supported Cell Types

Currently supported cell types (based on available datasets):
- **PBMC** - Peripheral Blood Mononuclear Cells
- **K562** - Chronic myelogenous leukemia cell line
- **HEK293** - Human Embryonic Kidney 293 cells
- **HCT116** - Human colorectal carcinoma cell line

## Data Sources

The script processes TFB datasets from three sources:
1. **UniBind** - Unified binding data
2. **ChIP-Atlas** - ChIP-seq data atlas
3. **ReMap2022** - Regulatory elements mapping

## Output

For each cell type and data source combination, the script generates:
- CSV files with GRN edges: `{CELL_TYPE}_{source}_{cell_type}.csv`
- Summary files: `tfb_{cell_type}_summary.txt`

Output files are saved to: `resources/grn_benchmark/ground_truth/`

## Examples

### Process K562 only:
```bash
cd /path/to/task_grn_inference
python src/process_data/tfb/process_tfb_datasets.py --cell-types k562
```

### Process PBMC and K562:
```bash
cd /path/to/task_grn_inference
python src/process_data/tfb/process_tfb_datasets.py --cell-types pbmc k562
```

### See what cell types are available:
```bash
cd /path/to/task_grn_inference
python src/process_data/tfb/process_tfb_datasets.py --list-available
```

## Requirements

- Python 3.x
- pandas
- Other dependencies from `build_grn.py`

The script expects to be run from the `task_grn_inference` root directory and will automatically locate the required resource files.