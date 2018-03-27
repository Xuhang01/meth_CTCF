import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np
import os
from concurrent import futures
######## probe position
meth_prot = '/f/hang/PhD_Project/Data/Database/ICGC/Methylation_Protocol_file/HumanMethylation450_15017482_v1-2.csv'
df_prot = pd.read_csv(meth_prot, skiprows=range(7))
df_prot = df_prot[df_prot.IlmnID.apply(lambda s: s.startswith('cg'))]
df_prot.loc[:,"order"] = df_prot.IlmnID.apply(lambda s: int(s[2:])/100000)
df_prot = df_prot.sort_values(by=["IlmnID"])
df_prot.loc[:,"CHR"] = df_prot.CHR.apply(lambda s: str(int(s)) if s not in ["X","Y","M"] else s )
df_prot.loc[:,"MAPINFO"] = df_prot.MAPINFO.astype(np.int)
df_prot.loc[:,["IlmnID","CHR","MAPINFO"]].sort_values(["CHR","MAPINFO"])\
.to_csv("./data/eQTL_input_icgc_donor_id/meth450_probe.txt",sep="\t",index=False)

######## gene position
ensembl_gene = "./data/ENSEMBL_gene.gtf"
col_dtypes = {
    "chr": 'category',
    "database": 'category',
    "type": 'category',
    "start": 'int',
    "end": 'int',
    "unknown": 'category',
    "strand": 'category',
    "unknown2": 'category',
    "info": 'object'
}
df_gene = pd.read_table(ensembl_gene, comment='#',
    names=["chr","database","type","start","end","unknown","strand","unknown2","info"],
    dtype=col_dtypes)
df_gene.loc[:,"gene_id"] = df_gene["info"].apply(lambda s: {each.split(' ')[0]: each.split(' ')[1] for each in s.split('; ')}["gene_name"].split('"')[1])
df_gene = df_gene[(df_gene.type=='CDS')]
df_gene = df_gene.groupby("gene_id").agg('first')
df_gene = df_gene.reset_index()

df_gene[["gene_id","chr","start","end"]].to_csv("./data/eQTL_input_icgc_donor_id/ensembl_gene.txt", sep="\t", index=False)

######## expression matrix
df_exp = pd.read_hdf("/f/hang/PhD_Project/Project/meth_anchor/data/exp_matrix/exp_matrix.ICGC.hdf",key="ICGC")
df_exp.columns = df_exp.columns.droplevel(level=0)
df_exp = df_exp.reset_index()
d_exp = df_exp.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
df_exp = df_exp.groupby("icgc_donor_id").agg(np.mean)


def generate_eQTL_input(df_exp, i_th):
    meth_matrix_dir = "/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/ICGC"
    df_meth = pd.read_hdf(os.path.join(meth_matrix_dir, "meth_matrix.ICGC_{}.hdf".format(i_th)), \
        key="ICGC_{}".format(i_th))
    df_meth = df_meth.reset_index()
    sample_type = "icgc_donor_id"
    df_meth = df_meth.drop(columns=["project_code","icgc_specimen_id","icgc_sample_id"])
    df_meth = df_meth.groupby("icgc_donor_id").agg(np.mean)
    df_meth = df_meth.dropna(axis=1)
    common_donor = set(df_exp.index.tolist()).intersection(df_meth.index.tolist())
    if len(common_donor) != 5455:
        print("ERROR")
    df_exp = df_exp.loc[common_donor, :]
    df_meth = df_meth.loc[common_donor, :]
    df_meth.T.to_csv("./data/eQTL_input_icgc_donor_id/meth_matrix_ith/meth_matrix_{}.txt".format(i_th),sep="\t")

    if os.path.exists("./data/eQTL_input_icgc_donor_id/exp_matrix.txt"):
        pass
    else:
        df_exp.T.to_csv("./data/eQTL_input_icgc_donor_id/exp_matrix.txt",sep="\t")

with futures.ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(generate_eQTL_input, [df_exp for each in range(277)], range(277))



# gene_list = df_exp.columns.tolist()
# probe_list = df_meth.columns.tolist()

# out_list = []
# for probe_id in probe_list[:10]:
#     for gene_id in gene_list[:10]:
#         sub_meth = df_meth.loc[:,[sample_type, probe_id]]
#         sub_gene = df_exp.loc[:,[sample_type, gene_id]]
#         sub_meth = sub_meth.groupby("icgc_donor_id").agg(np.mean)
#         sub_gene = sub_gene.groupby("icgc_donor_id").agg(np.mean)
#         df = pd.concat([sub_meth, sub_gene], axis=1)

#         fit = sm.OLS(df[gene_id], df[probe_id], missing='drop').fit()

#         a = fit.outlier_test()
#         df.loc[a.index,:]
#         a1 = a[a.iloc[:,2]<0.05]
#         df_new = df.loc[a1.index, :]
#         new_fit = sm.OLS(df_new[gene_id], df_new[probe_id], missing='drop').fit()
        
#         out_list.append([probe_id, gene_id, new_fit.pvalues])
#         if new_fit.pvalues[0] < 0.000001 :
#             print probe_id, gene_id, new_fit.pvalues
#             plt.scatter(df[probe_id], df[gene_id])
#             plt.xlim(0,1)
#             plt.ylim(0,max(sub_gene[gene_id]))
#             plt.show()

