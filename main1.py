from concurrent import futures
from models.DMS_model import *
from utils.common import *
import utils.common as common
print(dir(common))
if __name__ == '__main__':
    ##### call DMS for one case
    # proj_code = 'BRCA-US'
    # i_th = 0
    # df_meth = load_meth(proj_code, i_th)
    # ##### load expression file
    # df_exp = load_exp(proj_code)
    # ##### load sample file
    # df_sp = icgc_specimen(proj_code)
    # #####
    # df_DMS = DMS(df_meth, df_sp)
    # df_DMS.to_csv('./misc/DMS/stat/{}/DMS_stat.{}_{}.tsv'.format(proj_code, proj_code, i_th), sep='\t',index=False)
    
    ##### recurrently call DMS with multi-threading
    with futures.ProcessPoolExecutor(max_workers=20)  as executor:
        executor.map(call_DMS, proj_list)
