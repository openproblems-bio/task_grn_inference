# How to Add a New GRN Inference Method to geneRNBI

---

Step-by-step guide to integrating a new GRN inference method into the geneRNBI benchmark, including file structure, required format, and how to test with Viash.


## Overview

Every GRN inference method in geneRNBI lives under `src/methods/<method_name>/` and consists of exactly at least two required files:
- `config.vsh.yaml` — Viash component configuration
- `script.py` (or `script.R`) — the inference logic

---

## Step 1: Create the directory

```bash
mkdir -p src/methods/<method_name>
```

Example for Pearson correlation:
```bash
mkdir -p src/methods/pearson_corr
```

---

## Step 2: Write `script.py`

**IMPORTANT — always write files using bash, not Python open().**
Writing Python source code that contains docstrings (`"""`) inside a Python `open()` call causes syntax errors (quote conflicts). Always use a bash heredoc instead:

```bash
cat > /absolute/path/to/script.py << 'PYEOF'
<file content here>
PYEOF
```

The script MUST follow this exact structure. The `## VIASH START` / `## VIASH END` block is mandatory — Viash replaces it at runtime with the actual CLI parameters.

```python
import sys
import anndata as ad
import numpy as np
import pandas as pd

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'tf_all': 'resources_test/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/prediction.h5ad'
}
## VIASH END

# Load data
rna = ad.read_h5ad(par['rna'])
tf_all = np.loadtxt(par['tf_all'], dtype=str)

# ── YOUR METHOD LOGIC HERE ──────────────────────────────────────────────────
# Example: compute something and produce a net DataFrame with columns:
#   source (TF name), target (gene name), weight (float score)

# net = pd.DataFrame({'source': [...], 'target': [...], 'weight': [...]})
# ────────────────────────────────────────────────────────────────────────────

# Filter to TF-gene pairs only (important for evaluation — only top 50k edges scored)
net = net[net['source'].isin(tf_all)]

# Save output as AnnData
net['weight'] = net['weight'].astype(str)
output = ad.AnnData(
    X=None,
    uns={
        'method_id': '<method_name>',
        'dataset_id': rna.uns.get('dataset_id', 'unknown'),
        'prediction': net[['source', 'target', 'weight']]
    }
)
out_dir = os.path.dirname(par['prediction'])
if out_dir:
    os.makedirs(out_dir, exist_ok=True)
output.write_h5ad(par['prediction'])
print(f"Saved {len(net)} edges to {par['prediction']}")
```

**Critical rules for `script.py`:**
- Use `output.write_h5ad(par['prediction'])` — NOT `output.write()`
- Guard `makedirs` against empty string: `if out_dir: os.makedirs(out_dir, exist_ok=True)`
- All inference logic must run at module level (not inside `if __name__ == '__main__':`), so Viash can execute it directly

### Pearson correlation example (complete script.py):

```python
import sys
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'tf_all': 'resources_test/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/prediction.h5ad'
}
## VIASH END

rna = ad.read_h5ad(par['rna'])
tf_all = np.loadtxt(par['tf_all'], dtype=str)

# Get expression matrix (cells x genes)
X = rna.layers.get('lognorm', rna.X)
if hasattr(X, 'toarray'):
    X = X.toarray()
genes = np.array(rna.var_names)

# Find TF indices in the expression matrix
tf_mask = np.isin(genes, tf_all)
tf_indices = np.where(tf_mask)[0]
tf_names = genes[tf_indices]

# Compute Pearson correlation between each TF and all genes
rows = []
for i, tf_idx in enumerate(tf_indices):
    tf_expr = X[:, tf_idx]
    for j, gene in enumerate(genes):
        if genes[j] == tf_names[i]:
            continue  # skip self-loops
        corr = np.corrcoef(tf_expr, X[:, j])[0, 1]
        rows.append({'source': tf_names[i], 'target': gene, 'weight': abs(corr)})

net = pd.DataFrame(rows)
net = net.sort_values('weight', ascending=False).head(50000)

net['weight'] = net['weight'].astype(str)
output = ad.AnnData(
    X=None,
    uns={
        'method_id': 'pearson_corr',
        'dataset_id': rna.uns.get('dataset_id', 'unknown'),
        'prediction': net[['source', 'target', 'weight']]
    }
)
import os
out_dir = os.path.dirname(par['prediction'])
if out_dir:
    os.makedirs(out_dir, exist_ok=True)
output.write_h5ad(par['prediction'])
print(f"Saved {len(net)} edges to {par['prediction']}")
```

---

## Step 3: Write `config.vsh.yaml`

```yaml
__merge__: /src/api/comp_method.yaml

name: <method_name>
namespace: "grn_methods"
info:
  label: "<Pretty Method Name>"
  summary: "<One sentence description>"
  description: |
    <Longer description of your method>

resources:
  - type: python_script
    path: script.py

engines:
  - type: docker
    image: ghcr.io/openproblems-bio/base_python:1.0.4
    __merge__: /src/api/base_requirements.yaml
    setup:
      - type: python
        packages: [ numpy, pandas, scipy, anndata ]  # add your packages here

runners:
  - type: executable
  - type: nextflow
    directives:
      label: [midtime, midmem, midcpu]
```

### Pearson correlation example (complete config.vsh.yaml):

```yaml
__merge__: /src/api/comp_method.yaml

name: pearson_corr
namespace: "grn_methods"
info:
  label: "Pearson Correlation"
  summary: "Infers GRN using Pearson correlation between TF and target gene expression."
  description: |
    Computes Pearson correlation coefficients between all transcription factors
    and all target genes using log-normalized expression data.

resources:
  - type: python_script
    path: script.py

engines:
  - type: docker
    image: ghcr.io/openproblems-bio/base_python:1.0.4
    __merge__: /src/api/base_requirements.yaml
    setup:
      - type: python
        packages: [ numpy, pandas, scipy, anndata ]

runners:
  - type: executable
  - type: nextflow
    directives:
      label: [midtime, midmem, midcpu]
```

---

## Step 4: Download test data (first time only)

```bash
aws s3 sync s3://openproblems-data/resources_test/grn/grn_benchmark \
    resources_test/grn_benchmark --no-sign-request
```

---

## Step 5: Test with Viash

ALWAYS use `subprocess.run` to capture output and detect errors:

```python
import subprocess
result = subprocess.run(
    ["viash", "test", "src/methods/<method_name>/config.vsh.yaml"],
    capture_output=True, text=True
)
print(result.stdout)
print(result.stderr)
if result.returncode != 0:
    print("TEST FAILED — see error above")
else:
    print("TEST PASSED")
```

Or from the command line:
```bash
viash test src/methods/<method_name>/config.vsh.yaml
```

---


## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `Cannot connect to Docker daemon` | Docker not running | Start Docker Desktop (macOS/Windows) or `sudo systemctl start docker` (Linux), then verify with `docker info` |
| `__merge__ file not found` | Running from wrong directory | Always run viash commands from the repo root |
| `KeyError: dataset_id` | h5ad missing metadata | Use `rna.uns.get('dataset_id', 'unknown')` |
| `Empty prediction` | TF filter removed all rows | Check `tf_all` path and that gene names match `rna.var_names` |
