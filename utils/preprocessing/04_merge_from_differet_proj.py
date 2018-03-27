#!/usr/bin/env python

#### import module
from __future__ import division
import pandas as pd
import numpy as np
from concurrent import futures

proj_list = ['BLCA-US','BRCA-US','CESC-US','COAD-US','GBM-US','KIRC-US','KIRP-US','LAML-US',
             'LIHC-US','LUAD-US','LUSC-US','OV-US','PAAD-US','PRAD-US','READ-US','SKCM-US',
             'STAD-US','THCA-US','UCEC-US']

matrix_dir = "/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/"
merged_matrix_dir = "/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/ICGC/"
##### merge_matrix for each i_th  (totally 277)
def merge_matrix(i_th):
    df_matrix_list = [pd.read_hdf(matrix_dir+"/{}/meth_matrix.{}_{}.hdf".format(proj_code, proj_code, i_th), 
        key="{}_{}".format(proj_code, i_th)) for proj_code in proj_list]
    df_ICGC_ith = pd.concat(df_matrix_list)
    ICGC_file = merged_matrix_dir+"/meth_matrix.ICGC_{}.hdf".format(i_th)
    df_ICGC_ith.to_hdf(ICGC_file, mode='w', key="ICGC_{}".format(i_th))

with futures.ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(merge_matrix, np.arange(277))