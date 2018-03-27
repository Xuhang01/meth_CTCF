def mem_usage(pandas_obj):
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb)


for dtype in ['float','int','object']:
   selected_dtype = df.select_dtypes(include=[dtype])
   mean_usage_b = selected_dtype.memory_usage(deep=True).mean()
   mean_usage_mb = mean_usage_b / 1024 ** 2
   print("Average memory usage for {} columns: {:03.2f} MB".format(dtype,mean_usage_mb))


df_int = df.select_dtypes(include=['int'])
df_float = df.select_dtypes(include=['float'])
df_obj = df.select_dtypes(include=['object'])


converted_int = df_int.apply(pd.to_numeric,downcast='unsigned')
converted_float = df_float.apply(pd.to_numeric,downcast='float')

# not always change it
converted_obj = pd.DataFrame()
for col in df_obj.columns:
    num_unique_values = len(df_obj[col].unique())
    num_total_values = len(df_obj[col])
    if num_unique_values / num_total_values < 0.5:
        converted_obj.loc[:,col] = df_obj[col].astype('category')
    else:
        converted_obj.loc[:,col] = df_obj[col]

optimized_df = df.copy()
optimized_df[converted_int.columns] = converted_int
optimized_df[converted_float.columns] = converted_float
optimized_df[converted_obj.columns] = converted_obj


for dtype in ['float','int','object']:
   selected_dtype = df.select_dtypes(include=[dtype])
   mean_usage_b = selected_dtype.memory_usage(deep=True).mean()
   mean_usage_mb = mean_usage_b / 1024 ** 2
   print("Average memory usage for {} columns: {:03.2f} MB".format(dtype,mean_usage_mb))

for dtype in ['float','int','object']:
   selected_dtype = optimized_df.select_dtypes(include=[dtype])
   mean_usage_b = selected_dtype.memory_usage(deep=True).mean()
   mean_usage_mb = mean_usage_b / 1024 ** 2
   print("Average memory usage for {} columns: {:03.2f} MB".format(dtype,mean_usage_mb))


import pandas as pd
def mem_usage(pandas_obj):
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb)
column_types = {'icgc_donor_id': 'category',
                'project_code': 'category',
                'icgc_specimen_id': 'category',
                'icgc_sample_id': 'category',
                'submitted_sample_id': 'category',
                'analysis_id': 'category',
                'array_platform': 'category',
                'probe_id': 'category',
                'methylation_value': 'float32',
                'metric_used': 'category',
                'methylated_probe_intensity': 'category',
                'unmethylated_probe_intensity': 'category',
                'verification_status': 'category',
                'verification_platform': 'category',
                'fraction_wg_cpg_sites_covered': 'category',
                'conversion_rate': 'category',
                'experimental_protocol': 'category',
                'other_analysis_algorithm': 'category',
                'raw_data_repository': 'category',
                'raw_data_accession': 'category'
            }

read_and_optimized = pd.read_table('meth_array.LUAD-US.tsv.gz',compression='gzip' ,dtype=column_types)
print(mem_usage(read_and_optimized))

df = pd.read_table('test.tsv')
print(mem_usage(df))
