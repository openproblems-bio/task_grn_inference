
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys
import os



from regression_2.helper import main as main_reg2
from ws_distance.helper import main as main_ws_distance
from sem.helper import main as main_sem


def main(par):
    try:
        rr_reg2 = main_reg2(par)
    except Exception as e:
        print(f"Error in regression 2 metrics: {e}")
        rr_reg2 = pd.DataFrame()
    try:
        rr_sem = main_sem(par)
    except Exception as e:
        print(f"Error in sem metrics: {e}")
        rr_sem = pd.DataFrame()

    try:
        _, rr_ws = main_ws_distance(par)
    except Exception as e:
        print(f"Error in ws distance metrics: {e}")
        rr_ws = pd.DataFrame()

    
    

    rr_all = pd.concat([rr_reg2, rr_ws, rr_sem], axis=1)
    assert rr_all.shape[1] >1
    return rr_all
