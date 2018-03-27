from __future__ import division
import pandas as pd
import numpy as np
import os 

# ##### argparser
# import argparse
# parser = argparse.ArgumentParser(description="plot gene - meth_probe pairs")
# parser.add_argument('--i_th', help='batch id')
# parser.add_argument('--probe_id', help='meth probe id')
# parser.add_argument('--gene_id', help='ensemble gene id')
# args = parser.parse_args()
# i_th = args.i_th
# probe_id = args.probe_id
# gene_id = args.gene_id


##### load exp matrix
df_exp = pd.read_hdf("/f/hang/PhD_Project/Project/meth_anchor/data/exp_matrix/exp_matrix.ICGC.hdf",key="ICGC")
df_exp.columns = df_exp.columns.droplevel(level=0)
df_exp = df_exp.reset_index()
df_exp = df_exp.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
df_exp = df_exp.groupby("icgc_donor_id").agg(np.mean)
zero_exp_dict = {gene_id: np.sum(df_exp.loc[:,gene_id]==0) for gene_id in df_exp.columns.values}



def generate_meth_list(i_th):
    df_meth = pd.read_hdf("/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/ICGC/meth_matrix.ICGC_{}.hdf".format(i_th))
    df_meth = df_meth.reset_index()
    df_meth = df_meth.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
    df_meth = df_meth.groupby("icgc_donor_id").agg(np.mean)
    return df_meth
df_meth_list = []
with futures.ProcessPoolExecutor(max_workers=20) as executor:
    for each in executor.map(generate_meth_list, np.arange(277)):
        df_meth_list.append(each)

# def generate_gene_probe_pair_file(df_exp, probe_id, gene_id, i_th, beta, qvalue):
#     output_file = "/f/hang/PhD_Project/Project/meth_anchor/data/eQTL_input_icgc_donor_id/significant_eQTL_pairs/{}_{}.txt".format(gene_id, probe_id)
#     sub_exp = df_exp.loc[:, gene_id]
#     # if sub_exp.value_counts()[0] > 0:
#     #     print(probe_id, gene_id, 'fail1')
#     #     return False
#     df_meth = pd.read_hdf("/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/ICGC/meth_matrix.ICGC_{}.hdf".format(i_th))
#     df_meth = df_meth.reset_index()
#     df_meth = df_meth.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
#     df_meth = df_meth.groupby("icgc_donor_id").agg(np.mean)
#     sub_meth = df_meth.loc[:, probe_id]
#     if sub_meth.values.std() >0.1:
#         df_pairs = pd.concat([sub_exp, sub_meth],axis=1,join='inner')
#         df_pairs = df_pairs.reset_index()
#         df_pairs.to_csv(output_file, sep="\t")
#         os.system("Rscript ./misc/probe_gene_pairs.R {} {} {} {}".format(probe_id, gene_id, beta, qvalue))
#         print(probe_id, gene_id, 'success')
#         return True
#     else:
#         print(probe_id, gene_id, 'fail2')
#         return False



##### load all cis_eQTL results
df_list = []
for i_th in range(277):
    df = pd.read_table("./data/eQTL_input_icgc_donor_id/results/{}_cis_eQTL.txt".format(i_th))
    df.loc[:,"i_th"] = i_th
    df_list.append(df)
df_total = pd.concat(df_list)
df_total.loc[:,"zero_exp_count"] = df_total.gene.apply(lambda s: zero_exp_dict[s])

sig_df = df_total[ (df_total.beta.abs()>0.01) & (df_total.FDR < 1e-5) & (df_total.zero_exp_count <1)]


# ##### visualization 
for i in range(1000):
    probe_id, gene_id, beta, stat, pvalue, qvalue, i_th, zero_exp_count = sig_df.iloc[i,:]
    output_file = "/f/hang/PhD_Project/Project/meth_anchor/data/eQTL_input_icgc_donor_id/significant_eQTL_pairs/{}_{}.txt".format(gene_id, probe_id)
    sub_exp = df_exp.loc[:, gene_id]
    sub_meth = df_meth_list[i_th].loc[:, probe_id]

    if sub_meth.values.std() > 0.1:
        df_pairs = pd.concat([sub_exp, sub_meth],axis=1,join='inner')
        df_pairs = df_pairs.reset_index()
        df_pairs.to_csv(output_file, sep="\t")
        os.system("Rscript ./misc/probe_gene_pairs.R {} {} {} {}".format(probe_id, gene_id, beta, qvalue))
        print(probe_id, gene_id, 'success')
    else:
        print(probe_id, gene_id, 'fail2')
