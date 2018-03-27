### ICGC exp_seq.proj_code.tsv.gz -> exp_matrix.proj_code.hdf
#   exp_matrix.proj_code.hdf:
#              column:   gene_id
#              multiindex:   [project_code, icgc_donor_id, icgc_specimen_id, icgc_sample_id]
import pandas as pd
import numpy as np
from concurrent import futures
column_types = {'icgc_donor_id': 'category',
                'project_code': 'category',
                'icgc_specimen_id': 'category',
                'icgc_sample_id': 'category',
                'submitted_sample_id': 'category',
                'analysis_id': 'category',
                'gene_model': 'category',
                'gene_id': 'category',
                'normalized_read_count': 'float32',
                'raw_read_count': 'float32',
                'fold_change': 'float32',
                'assembly_version': 'category',
                'platform': 'category',
                'total_read_count': 'float32',
                'experimental_protocol': 'category',
                'alignment_algorithm': 'category',
                'normalization_algorithm': 'category',
                'other_analysis_algorithm': 'category',
                'sequencing_strategy': 'category',
                'raw_data_repository': 'category',
                'raw_data_accession': 'category',
                'reference_sample_type': 'category'
            }

proj_list = ['BLCA-US','BRCA-US','CESC-US','COAD-US','GBM-US','KIRC-US','KIRP-US','LAML-US',
                     'LIHC-US','LUAD-US','LUSC-US','OV-US','PAAD-US','PRAD-US','READ-US','SKCM-US','STAD-US','THCA-US','UCEC-US']

##### transfer expression dataframe to matrix for each project
def write_exp_matrix(proj):
    exp_file = '/f/hang/PhD_Project/Data/Database/ICGC/ICGC_Project/{}/exp_seq.{}.tsv.gz'.format(proj,proj)
    df_exp = pd.read_table(exp_file, compression='gzip', dtype=column_types)
    ## raw_read_count to nomralized_reads_count
    # raw_read_count/third_quantile*1000  (3rd quantile in each sample)
    new_list = []
    for sample_id in df_exp.icgc_sample_id.values.unique():
        sub_exp = df_exp[df_exp.icgc_sample_id == sample_id]
        quantile = sub_exp[sub_exp.raw_read_count != 0].raw_read_count.quantile(0.75)
        sub_exp.loc[:,'quantile'] = quantile
        sub_exp.loc[:,"normalized_count"] = (sub_exp.raw_read_count+1)/quantile*1000
        sub_exp.loc[:,"log2_normalized_count"] = np.log2(sub_exp.normalized_count)
        new_list.append(sub_exp)
    ## concatenate expression files
    new_df_exp = pd.concat(new_list)
    ## select some needed columns
    new_df_exp = new_df_exp.loc[:,["project_code","icgc_donor_id","icgc_specimen_id","icgc_sample_id","gene_id","log2_normalized_count"]]
    new_df_exp = new_df_exp.drop_duplicates(subset=["project_code","icgc_donor_id","icgc_specimen_id","icgc_sample_id","gene_id"])
    ## transfer dataframe to matrix: columns, index
    exp_matrix = new_df_exp.set_index(["project_code","icgc_donor_id","icgc_specimen_id","icgc_sample_id","gene_id"])
    exp_matrix = exp_matrix.unstack('gene_id')
    exp_matrix.columns = exp_matrix.columns.droplevel(level=0)
    exp_matrix_file = "./data/exp_matrix/exp_matrix.{}.hdf".format(proj)
    exp_matrix.to_hdf(exp_matrix_file, mode='w', key=proj)

with futures.ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(write_exp_matrix, proj_list)


###### merge all exp matrix files
df_ICGC_matrix = pd.concat([
    pd.read_hdf("./data/exp_matrix/"+"exp_matrix.{}.hdf".format(proj))
    for proj in proj_list
    ])
df_ICGC_matrix.to_hdf("./data/exp_matrix/exp_matrix.ICGC.hdf", 
    mode='w',key="ICGC")