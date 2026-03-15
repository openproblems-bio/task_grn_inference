# Viash Commands for GeneRNBI

---

## Metadata

**Short Description**: Reference guide for all Viash CLI commands used in GeneRNBI — running, testing, and building GRN inference methods and evaluation metrics.

**Authors**: GeneRNBI Team

**Version**: 1.0

**Last Updated**: March 2026

**License**: CC BY 4.0

**Commercial Use**: ✅ Allowed

---

## Overview

Viash is the component framework used by GeneRNBI to package and run GRN methods. All Viash commands must be run from the **repository root** (`task_grn_inference/`).

CRITICAL: Always capture command output using `subprocess.run(capture_output=True, text=True)` so errors are visible.

---

## Running a Method

### Run with default test parameters (executable runner, no Docker):
```bash
viash run src/methods/<method_name>/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --prediction output/net.h5ad \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv
```

### Run a specific dataset:
```bash
viash run src/methods/pearson_corr/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/norman_rna.h5ad \
    --prediction output/norman_net.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv
```

---

## Testing a Method (requires Docker)

### Test (builds Docker image and runs):
```bash
viash test src/methods/<method_name>/config.vsh.yaml
```

### Python equivalent (always use this to capture errors):
```python
import subprocess
result = subprocess.run(
    ["viash", "test", "src/methods/pearson_corr/config.vsh.yaml"],
    capture_output=True, text=True, cwd="/path/to/repo/root"
)
print(result.stdout)
print(result.stderr)
if result.returncode != 0:
    print("FAILED")
else:
    print("PASSED")
```

---

## Building a Docker Image

```bash
viash build src/methods/<method_name>/config.vsh.yaml -o output/ --engine docker
```

---

## Inspecting a Config

```bash
viash config view src/methods/<method_name>/config.vsh.yaml
```

---

## Running Without Docker (Executable Runner)

When Docker is unavailable, use the executable runner directly:
```bash
# First build the executable
viash build src/methods/<method_name>/config.vsh.yaml -o output/ --engine native

# Then run it
output/<method_name> --rna <rna_file> --prediction <output_file> --tf_all <tf_file>
```

---

## Running All Methods (Pipeline)

```bash
nextflow run main.nf -params-file params/params.yaml
```

---

## Troubleshooting

### `Cannot connect to Docker daemon`
Docker is not running. Options:
1. Start Docker: `systemctl start docker` (or Docker Desktop)
2. Use native/executable runner instead of test: `viash run` instead of `viash test`
3. Check socket: `export DOCKER_HOST=unix://$XDG_RUNTIME_DIR/docker.sock`

### `__merge__ file not found: /src/api/comp_method.yaml`
You are not running from the repo root. Always `cd` to `task_grn_inference/` first.

### `Error: Unknown argument`
The `--` separator between viash args and component args is required:
```bash
viash run config.vsh.yaml -- --rna file.h5ad   # correct
viash run config.vsh.yaml --rna file.h5ad       # wrong
```

### `Setup failed` during `viash test`
Docker build failed. Check:
- Package names are correct in `config.vsh.yaml` setup section
- Base image exists and is accessible
- Internet access for pip install

---

## Config.vsh.yaml Format Quick Reference

```yaml
__merge__: /src/api/comp_method.yaml    # REQUIRED — merges common method schema

name: my_method                          # unique snake_case ID
namespace: "grn_methods"                 # always "grn_methods" for inference methods

info:
  label: "My Method"                     # human-readable name
  summary: "Short description"

resources:
  - type: python_script                  # or r_script
    path: script.py

engines:
  - type: docker
    image: ghcr.io/openproblems-bio/base_python:1.0.4
    __merge__: /src/api/base_requirements.yaml
    setup:
      - type: python
        packages: [ numpy, pandas ]      # extra pip packages

runners:
  - type: executable
  - type: nextflow
    directives:
      label: [midtime, midmem, midcpu]   # resource label
```

Resource labels (from `scripts/labels_tw.config`):
- `[lowtime, lowmem, lowcpu]` — quick/small jobs
- `[midtime, midmem, midcpu]` — typical methods
- `[hightime, highmem, highcpu]` — heavy methods (e.g. deep learning)
