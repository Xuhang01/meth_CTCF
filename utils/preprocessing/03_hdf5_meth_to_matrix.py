import pandas as pd
import numpy as np
from concurrent import futures
columns = ['icgc_donor_id','project_code','icgc_specimen_id','icgc_sample_id',
           'submitted_sample_id','analysis_id','array_platform','probe_id',
           'methylation_value','metric_used','methylated_probe_intensity',
           'unmethylated_probe_intensity','verification_status','verification_platform',
           'fraction_wg_cpg_sites_covered','conversion_rate','experimental_protocol',
           'other_analysis_algorithm','raw_data_repository','raw_data_accession'
        ]
proj_list = ['BLCA-US','BRCA-US','CESC-US','COAD-US','GBM-US','KIRC-US','KIRP-US','LAML-US',
             'LIHC-US','LUAD-US','LUSC-US','OV-US','PAAD-US','PRAD-US','READ-US','SKCM-US',
             'STAD-US','THCA-US','UCEC-US']
##### process protocol file
meth_prot = '/f/hang/PhD_Project/Data/Database/ICGC/Methylation_Protocol_file/HumanMethylation450_15017482_v1-2.csv'
df_prot = pd.read_csv(meth_prot, skiprows=range(7))
df_prot = df_prot[df_prot.IlmnID.apply(lambda s: s.startswith('cg'))]
df_prot.loc[:,"order"] = df_prot.IlmnID.apply(lambda s: int(s[2:])//100000)
df_prot = df_prot.sort_values(by=["IlmnID"])
# split protocol
probe_split = []
for i in range(277):
    probe_split.append(df_prot[df_prot.order==i].IlmnID.tolist())
##### process protocol file
def raw_to_matrix(proj_code, i_th):
    #### raw_file in ./data/split_raw_new
    raw_file_dir = "./data/split_raw/"
    raw_meth = raw_file_dir + "/{}/{}_{}.tsv.gz".format(proj_code, proj_code, i_th)
    #### read raw file
    # df_meth = pd.read_table(raw_meth, compression='gzip',dtype=column_types)
    df_meth = pd.read_table(raw_meth, compression='gzip', dtype=column_types, names=columns)
    #### create matrix with either pivot_table or set_index followed by unstack
    ## pivot_table: to slow because the usage of four level of index
    # df_matrix = pd.pivot_table(df_meth, values="methylation_value", columns="probe_id", 
    #     index=["project_code","icgc_donor_id","icgc_specimen_id","icgc_sample_id"]).dropna(axis=0,how='all')
    ## set_index and unstack, much faster
    df_meth = df_meth.loc[:,["project_code","icgc_donor_id","icgc_specimen_id","icgc_sample_id","probe_id","methylation_value"]]
    df_matrix = df_meth.set_index(["project_code","icgc_donor_id","icgc_specimen_id","icgc_sample_id","probe_id"]).unstack("probe_id")
    df_matrix.columns = df_matrix.columns.droplevel()
    #### fill the blank according to protocol library. Some probe maybe not included in the data
    df_matrix.columns = df_matrix.columns.astype('object')
    difference_probe = set(probe_split[i_th]).difference(set(df_matrix.columns.tolist()))
    for each in difference_probe:
        df_matrix.loc[:,each] = np.nan
    df_matrix = df_matrix.T.sort_index().T
    #### write to hdf file
    matrix_split_file = "./data/meth_matrix_hdf/{}/meth_matrix.{}_{}.hdf".format(proj_code,proj_code, i_th)
    df_matrix.to_hdf(matrix_split_file, mode='w', key="{}_{}".format(proj_code, i_th)) 

for i in range(277):
    with futures.ProcessPoolExecutor(max_workers=20) as executor:
        executor.map(raw_to_matrix, proj_list, [i for j in range(len(proj_list))])
