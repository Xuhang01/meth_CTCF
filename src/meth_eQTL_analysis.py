from __future__ import division
import pandas as pd
import numpy as np
import os 

import argparse

parser = argparse.ArgumentParser(description="plot gene - meth_probe pairs")
parser.add_argument('--i_th', help='batch id')
parser.add_argument('--probe_id', help='meth probe id')
parser.add_argument('--gene_id', help='ensemble gene id')

args = parser.parse_args()
i_th = args.i_th
probe_id = args.probe_id
gene_id = args.gene_id


eQTL_dir = "/f/hang/PhD_Project/Project/meth_anchor/data/eQTL_input_icgc_donor_id/results/"
i_th = 0
cis_eQTL_file =  os.path.join(eQTL_dir, "{}_cis_eQTL.txt".format(i_th))
df_cis_eQTL = pd.read_table(cis_eQTL_file)

def generate_gene_probe_pair_file(gene_id, probe_id, i_th):
    output_file = "/f/hang/PhD_Project/Project/meth_anchor/data/eQTL_input_icgc_donor_id/significant_eQTL_pairs/{}_{}.txt".format(gene_id, probe_id)
    df_exp = pd.read_hdf("/f/hang/PhD_Project/Project/meth_anchor/data/exp_matrix/exp_matrix.ICGC.hdf",key="ICGC")
    df_exp.columns = df_exp.columns.droplevel(level=0)
    df_exp = df_exp.reset_index()
    df_exp = df_exp.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
    df_exp = df_exp.groupby("icgc_donor_id").agg(np.mean)

    sub_exp = df_exp.loc[:, gene_id]

    df_meth = pd.read_hdf("/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/ICGC/meth_matrix.ICGC_{}.hdf".format(i_th))
    df_meth = df_meth.reset_index()
    df_meth = df_meth.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
    df_meth = df_meth.groupby("icgc_donor_id").agg(np.mean)

    sub_meth = df_meth.loc[:, probe_id]

    df_pairs = pd.concat([sub_exp, sub_meth],axis=1,join='inner')
    df_pairs = df_pairs.reset_index()
    df_pairs.to_csv(output_file, sep="\t")

# gene_id = "APOA2"
# i_th = 0
# probe_id = "cg00058923"

generate_gene_probe_pair_file(gene_id, probe_id, i_th)

os.system("Rscript ./misc/probe_gene_pairs.R {} {}".format(gene_id, i_th, probe_id))


