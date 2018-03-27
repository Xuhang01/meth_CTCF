from __future__ import division
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
######### constants
proj_list = ['BLCA-US','BRCA-US','CESC-US','COAD-US','GBM-US','KIRC-US','KIRP-US','LAML-US',
             'LIHC-US','LUAD-US','LUSC-US','OV-US','PAAD-US','PRAD-US','READ-US','SKCM-US',
             'STAD-US','THCA-US','UCEC-US']
######### load data
meth_matrix_dir = "/f/hang/PhD_Project/Project/meth_anchor/data/meth_matrix_hdf/ICGC/"
i_th = 0
meth_matrix_ith = os.path.join(meth_matrix_dir, "meth_matrix.ICGC_{}.hdf".format(i_th))
df_matrix = pd.read_hdf(meth_matrix_ith, key="ICGC_{}".format(i_th))


######## 1. generate donor_count.txt
# count donor, specimen, sample
donor_count = []
specimen_count = []
sample_count = []
for proj_code in proj_list:
    sub_matrix = df_matrix[df_matrix.index.get_level_values('project_code') == proj_code]
    donor_count.append(len(sub_matrix.index.get_level_values('icgc_donor_id').unique()))
    specimen_count.append(len(sub_matrix.index.get_level_values('icgc_specimen_id').unique()))
    sample_count.append(len(sub_matrix.index.get_level_values('icgc_sample_id').unique()))
df_count = pd.DataFrame(
    data = np.array([donor_count, specimen_count, sample_count]), 
    columns = proj_list,
    index = ["icgc_donor_id","icgc_specimen_id","icgc_sample_id"]
    ).T
df_count.to_csv("./misc/stat/donor_count.txt", sep='\t')

######## 2. generate separate files for each probe_id (of course not so many)
def write_to_separate(df_matrix, probe_id):
    sub_matrix = df_matrix.loc[:,probe_id]
    sub_matrix = sub_matrix.reset_index()
    sub_matrix.to_csv("./misc/stat/probe_id/{}.txt".format(probe_id),sep='\t',index=False)

probe_id_list = df_matrix.columns.tolist()[:10]
for probe_id in probe_id_list:
    write_to_separate(df_matrix, probe_id)

######## 3. mean and std
mean_array = df_matrix.mean(axis=0).values
std_array = df_matrix.std(axis=0).values
df_stat = pd.DataFrame(data=np.array([mean_array, std_array]).T, columns=["mean_meth","std_meth"],index=df_matrix.columns)
df_stat.to_csv("./misc/stat/mean_and_std.txt", sep='\t')


######## 4. Island_df and non_Island_df
import pandas as pd
import numpy as np
#### process protocol file
meth_prot = '/f/hang/PhD_Project/Data/Database/ICGC/Methylation_Protocol_file/HumanMethylation450_15017482_v1-2.csv'
df_prot = pd.read_csv(meth_prot, skiprows=range(7),dtypes=)
df_prot = df_prot[df_prot.IlmnID.apply(lambda s: s.startswith('cg'))]
df_prot.loc[:,"order"] = df_prot.IlmnID.apply(lambda s: int(s[2:])/100000)
df_prot = df_prot.sort_values(by=["IlmnID"])

#### probe related to Island or not
sub_df = df_prot[df_prot.order==0]
sub_df = sub_df[["IlmnID","CHR","MAPINFO"]]
Island_df = sub_df[~ ((pd.isnull(sub_df.UCSC_CpG_Islands_Name)) & pd.isnull(sub_df.HMM_Island))]
non_Island_df = sub_df[ ((pd.isnull(sub_df.UCSC_CpG_Islands_Name)) & pd.isnull(sub_df.HMM_Island))]

#### check empty (no data at all in all sample)
case1 = 0
case2 = 0
for probe_id in df_matrix.columns.tolist():
    if probe_id in Island_df.IlmnID.tolist():
        if df_matrix[probe_id].isnull().sum() == df_matrix.shape[0]:
            case1+=1
        else:
            pass
    else:
        if df_matrix[probe_id].isnull().sum() == df_matrix.shape[0]:
            case2 +=1
        else:
            pass
## three comparison between Island and non-Island: case_count, mean_meth, std_meth 
case_count = [df_matrix[probe_id].notnull().sum() for probe_id in df_matrix.columns]
mean_meth = [df_matrix[probe_id].mean() for probe_id in df_matrix.columns]
std_meth = [df_matrix[probe_id].std() for probe_id in df_matrix.columns]
Island_or_not = [ 'Island' if probe_id in Island_df.IlmnID.tolist() else 'non_Island' for probe_id in df_matrix.columns]
df_out = pd.DataFrame(data=np.array([df_matrix.columns.tolist(), case_count, mean_meth, std_meth, Island_or_not]).T, columns=["probe_id","case_count","mean_meth","std_meth","is_Island"])
df_out.to_csv("./misc/stat/Island_and_non_Island_comparison.txt",sep='\t',index=False)

######## 5. gene correlation
