##### model to detect differentially methylated sites
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests
from utils.common import *
from concurrent import futures

##### call DMS with t_test
# NT: Nomral Tissue
# PTS: Primary Tumor Solid Tissue
# for each probe, 20% low NT vs 20% low PTS, 20% high NT vs 20% high PTS, using t-test to call DMS
# input: 
# output: 
def DMS(df_meth, df_sp):
    proj_code = df_sp.project_code.iat[0]
    # get specimen_type_dict
    specimen_list = df_sp.icgc_specimen_id.values
    type_list = df_sp.type.values
    specimen_type_dict = dict(zip(specimen_list, type_list))
    # add type to multiindex of specimen_type_dict
    df_meth.loc[:,"type"] = df_meth.index.get_level_values('icgc_specimen_id').to_series().apply(lambda s: specimen_type_dict[s]).tolist()
    df_meth.set_index('type',append=True,inplace=True)
    df_meth = df_meth.dropna(axis=1,how='all')
    meth_N = df_meth[df_meth.index.get_level_values('type').isin(['NB','NT'])]
    meth_T = df_meth[df_meth.index.get_level_values('type').isin(['PTS','MT'])]
    # for each probe_id, compare the 20% low and 20% high with ttest_ind (scipy.stats)
    # result = ttest_ind(meth_N, meth_T, axis=0, equal_var = False, nan_policy = 'omit')
    result_low = []
    result_high = []
    for probe_id in df_meth.columns.values:
        probe_N = meth_N[probe_id]
        probe_T = meth_T[probe_id]
        ### remove NaN
        probe_N = probe_N.dropna()
        probe_T = probe_T.dropna()
        ### low
        low_20_N = probe_N[probe_N.argsort().argsort()<len(probe_N)*0.2]
        low_20_T = probe_T[probe_T.argsort().argsort()<len(probe_T)*0.2]
        t_stats, pvalue = ttest_ind(low_20_T, low_20_N, axis=0, equal_var=False, nan_policy='omit')
        # if type(t_stats) == np.ma.core.MaskedConstant:
        #     t_stats = np.nan
        #     pvalue = np.nan
        effect_size = low_20_T.mean()-low_20_N.mean()
        result_low.append([probe_id, t_stats, pvalue, effect_size])
        ### high
        high_20_N = probe_N[probe_N.argsort().argsort()>len(probe_N)*0.8]
        high_20_T = probe_T[probe_T.argsort().argsort()>len(probe_T)*0.8]
        t_stats, pvalue = ttest_ind(high_20_N, high_20_T, axis=0, equal_var=False, nan_policy='omit')
        # if type(t_stats) == np.ma.core.MaskedConstant:
        #     t_stats = np.nan
        #     pvalue = np.nan
        effect_size = high_20_T.mean()-high_20_N.mean()
        result_high.append([probe_id, t_stats, pvalue, effect_size])
    ### low
    df_low = pd.DataFrame(result_low, columns=['probe_id','t_stats','pvalue','effect_size'],index=df_meth.columns)
    df_low = df_low.dropna(subset=['pvalue'],axis=0)
    p_reject, p_adjust, a, b = multipletests(df_low.pvalue.values, alpha=0.01, method='fdr_bh', is_sorted=False, returnsorted=False)
    df_low.loc[:,'p_reject'] = p_reject
    df_low.loc[:,'p_adjust'] = p_adjust
    ### high
    df_high = pd.DataFrame(result_high, columns=['probe_id','t_stats','pvalue','effect_size'],index=df_meth.columns)
    df_high = df_high.dropna(subset=['pvalue'],axis=0)
    p_reject, p_adjust, a, b = multipletests(df_high.pvalue.values, alpha=0.01, method='fdr_bh', is_sorted=False, returnsorted=False)
    df_high.loc[:,'p_reject'] = p_reject
    df_high.loc[:,'p_adjust'] = p_adjust
    ### return DMS
    df_low.loc[:,'type'] = 'low'
    df_high.loc[:, 'type'] = 'high'
    df_DMS = pd.concat([df_low, df_high])
    df_DMS.loc[:,'proj_code'] = proj_code
    return df_DMS

##### call DMS for one case
# proj_code = 'BRCA-US'
# i_th = 0
# df_meth = load_meth(proj_code, i_th)
# ##### load expression file
# df_exp = load_exp(proj_code)
# ##### load sample file
# df_sp = icgc_specimen(proj_code)
# #####
# df_DMS = DMS(df_meth, df_sp, i_th)
# df_DMS.to_csv('./misc/DMS/stat/{}/DMS_stat.{}_{}.tsv'.format(proj_code, proj_code, i_th), sep='\t',index=False)


##### recurrently call DMS
# for proj_code in proj_list:
#     df_exp = load_exp(proj_code)
#     df_sp = icgc_specimen(proj_code)
#     for i_th in range(277):
#         df_meth = load_meth(proj_code, i_th)
#         df_DMS = DMS(df_meth, df_sp)
#         df_DMS.to_csv('./misc/DMS/stat/{}/DMS_stat.{}_{}.tsv'.format(proj_code, proj_code, i_th), sep='\t',index=False)

##### recurrently call DMS with multi-threading
def call_DMS(proj_code):
    df_exp = load_exp(proj_code)
    df_sp = icgc_specimen(proj_code)
    for i_th in range(277):
        df_meth = load_meth(proj_code, i_th)
        df_DMS = DMS(df_meth, df_sp)
        df_DMS.to_csv('./misc/DMS/stat/{}/DMS_stat.{}_{}.tsv'.format(proj_code, proj_code, i_th), sep='\t',index=False)

if __name__ == '__main__':
    with futures.ProcessPoolExecutor(max_workers=20)  as executor:
        executor.map(call_DMS, proj_list)