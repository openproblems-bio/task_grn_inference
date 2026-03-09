**Short Description**: Common errors and troubleshooting guidance for GeneRNBI users.

---

## Docker daemon not running

**Error message:**
```
Cannot connect to the Docker daemon at unix:///.../.docker/run/docker.sock. Is the docker daemon running?
ERROR! Setup failed!
```

**Cause:** Viash uses Docker by default to build and run components. If the Docker daemon is not running, all `viash run`, `viash test`, and `viash build` commands that use the default `docker` engine will fail.

**Fix — start Docker:**
- macOS/Windows: open Docker Desktop and wait for it to finish starting.
- Linux: run `sudo systemctl start docker` (or `sudo service docker start`).
- Verify it is running: `docker info`

**Alternative — run without Docker:**
Docker is not required for all parts of the framework. If Docker is unavailable or not supported on your system, you can still:

1. **Run GRN evaluation** using the local script (bypasses Docker entirely):
   ```bash
   bash src/metrics/all_metrics/run_local.sh \
       --dataset <dataset_name> \
       --prediction=<inferred_GRN.h5ad> \
       --score <output_score_file.h5ad> \
       --num_workers <n>
   ```
   Example:
   ```bash
   bash src/metrics/all_metrics/run_local.sh \
       --dataset op \
       --prediction=resources/grn_models/op/collectri.h5ad \
       --score=output_score_file.h5ad \
       --num_workers=20
   ```

2. **Run GRN inference** directly in Python without Viash by calling the `infer_grn()` function from the method's `script.py` with the appropriate input paths.

3. **Evaluate a new GRN model** not in GeneRNBI — generate the consensus prior first:
   ```bash
   bash scripts/prior/run_consensus.sh --dataset op --new_model {new_grn_model_file.h5ad}
   ```

**See also:** `docs/source/evaluation.rst` → "Running GRN evaluation without docker" for the full guide.
