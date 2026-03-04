# Running and Testing Methods with Viash

---

## Metadata

**Short Description**: How to use Viash CLI to test, run, and build GeneRNBI methods — including correct subprocess usage patterns for capturing output and detecting errors.

**Authors**: GeneRNBI Team

**Version**: 1.0

**Last Updated**: March 2026

**License**: CC BY 4.0

**Commercial Use**: ✅ Allowed

---

## Overview

Viash wraps your `script.py` in a Docker container and wires up CLI parameters from `config.vsh.yaml`. The three most common commands are `viash test`, `viash run`, and `viash build`.

**Always use `subprocess.run(capture_output=True)` — never `os.system()`.** The latter discards stdout/stderr so errors are invisible to you.

---

## Correct Pattern for Running Viash Commands

```python
import subprocess

def run_viash(args: list, cwd: str = None) -> bool:
    """Run a viash command and print output. Returns True if successful."""
    result = subprocess.run(
        args,
        capture_output=True,
        text=True,
        cwd=cwd
    )
    print(result.stdout)
    print(result.stderr)
    if result.returncode != 0:
        print(f"❌ FAILED (exit code {result.returncode})")
        return False
    print("✅ Success")
    return True
```

---

## Common Commands

### Test a method (builds Docker image + runs unit tests)
```python
run_viash(["viash", "test", "src/methods/my_method/config.vsh.yaml"])
```
Requires Docker daemon running. If Docker is not available, use `viash run` instead.

### Run a method directly (no Docker build, uses executable runner)
```python
run_viash([
    "viash", "run", "src/methods/my_method/config.vsh.yaml", "--",
    "--rna", "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
    "--prediction", "output/net.h5ad",
    "--tf_all", "resources_test/grn_benchmark/prior/tf_all.csv"
])
```

### Build a Docker image
```python
run_viash(["viash", "build", "src/methods/my_method/config.vsh.yaml", "-o", "target/"])
```

---

## Detecting Docker Issues

If `viash test` output contains any of these, Docker is not running:
- `Cannot connect to the Docker daemon`
- `Setup failed!`
- `Is the docker daemon running?`

In that case, switch to `viash run` with the executable runner (no Docker needed):

```python
result = subprocess.run(
    ["viash", "run", "src/methods/my_method/config.vsh.yaml", "--",
     "--rna", "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
     "--prediction", "output/net.h5ad",
     "--tf_all", "resources_test/grn_benchmark/prior/tf_all.csv"],
    capture_output=True, text=True
)
print(result.stdout)
print(result.stderr)
```

---

## Viash Version and Location

```bash
viash --version   # should print 0.9.x
which viash       # /usr/local/bin/viash or similar
```

The GeneRNBI repo already has Viash configured. Do not install a new version.

---

## Writing Method Files to the Correct Location

When creating a new method, always write files to `src/methods/<method_name>/`, NOT to the current working directory:

```python
import os

method_name = "pearson_corr"
method_dir = f"src/methods/{method_name}"
os.makedirs(method_dir, exist_ok=True)

with open(f"{method_dir}/script.py", "w") as f:
    f.write(script_content)

with open(f"{method_dir}/config.vsh.yaml", "w") as f:
    f.write(config_content)
```
