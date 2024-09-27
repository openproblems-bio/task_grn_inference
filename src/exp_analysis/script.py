import pandas as pd
import anndata as ad
import sys
import json

## VIASH START
par = {
}
## VIASH END
meta = {
  "resources_dir":'src/exp_analysis/'
}
sys.path.append(meta["resources_dir"])

from main import main 
main(par)



