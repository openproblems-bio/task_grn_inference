"""Build small input files (op, norman) for process_datasets from existing grn_benchmark files."""
import shutil, sys
from pathlib import Path
import anndata as ad, numpy as np

FULL, TEST, OUT = Path('resources/grn_benchmark'), Path('resources_test/grn_benchmark'), Path(sys.argv[1] if len(sys.argv) > 1 else 'resources_test/task_grn_inference/inputs')
N_CELLS = 2000
rng = np.random.default_rng(0)


def subset(src, dst, n=N_CELLS):
    a = ad.read_h5ad(src)
    if a.n_obs > n:
        a = a[np.sort(rng.choice(a.n_obs, n, replace=False))].copy()
    dst.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(dst, compression='gzip')


def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)


# op: the existing test files are already small
for k, p in dict(rna='inference_data/op_rna', atac='inference_data/op_atac', evaluation='evaluation_data/op_bulk').items():
    copy(TEST / f'{p}.h5ad', OUT / 'op' / f'{k}.h5ad')

# norman: subset the full files
subset(FULL / 'inference_data/norman_rna.h5ad', OUT / 'norman/rna.h5ad')
copy(FULL / 'evaluation_data/norman_bulk.h5ad', OUT / 'norman/evaluation.h5ad')
copy(FULL / 'evaluation_data/norman_de.h5ad', OUT / 'norman/evaluation_de.h5ad')
subset(Path('resources/processed_data/norman_evaluation_sc.h5ad'), OUT / 'norman/evaluation_sc.h5ad')
copy(TEST / 'prior/ws_distance_background_norman.csv', OUT / 'norman/ws_background.csv')
