
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
from vc_v2.helper import main as main_vc_v2
from tf_recovery import main as main_tf_rec


def main(par):
    rr_store = []

    if True:
        try:
            tf_rec = main_tf_rec(par)
        except Exception as e:
            print(f"Error in main_tf_rec metrics: {e}")
            tf_rec = pd.DataFrame()
        print("tf_rec done: ", tf_rec)
        rr_store.append(tf_rec)

    
    # if True:
    #     try:
    #         rr_vc = main_vc_v2(par)
    #     except Exception as e:
    #         print(f"Error in vc metrics: {e}")
    #         rr_vc = pd.DataFrame()
    #     print("vc done: ", rr_vc)
    #     rr_store.append(rr_vc)


    # try:
    #     rr_reg2 = main_reg2(par)
    # except Exception as e:
    #     print(f"Error in regression 2 metrics: {e}")
    #     rr_reg2 = pd.DataFrame()
    # rr_store.append(rr_reg2)
    # print("reg2 done: ", rr_reg2)

    # if False:
    #     try:
    #         rr_sem = main_sem(par)
    #     except Exception as e:
    #         print(f"Error in sem metrics: {e}")
    #         rr_sem = pd.DataFrame()
    #     print("sem done: ", rr_sem)
    #     rr_store.append(rr_sem)
    
    
    # try:
    #     _, rr_ws = main_ws_distance(par)
    # except Exception as e:
    #     print(f"Error in ws distance metrics: {e}")
    #     rr_ws = pd.DataFrame()
    # print("ws done: ", rr_ws)
    # rr_store.append(rr_ws)

    
    

    rr_all = pd.concat(rr_store, axis=1)
    assert rr_all.shape[1] >0
    return rr_all
