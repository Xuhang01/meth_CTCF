import pandas as pd
import numpy as np

import scipy.stats.mannwhitneyu as mwu

def eQTM(df_probe, df_gene, test='mwu'):
    # df_probe: methylation of one probe and one project
    # df_gene:  exp of one gene and one project
    # df_probe:   icgc_specimen_id, meth_value, Tumor|Normal
    # df_gene:    icgc_specimen_id, exp_value, Tumor|Normal
    df_probe = df_probe.reset_index().iloc[:,[2,4]]
    df_probe.columns = ['icgc_specimen_id','meth_value']
    df_gene = df_gene.reset_index().iloc[:,[2,4]]
    df_gene.columns = ['icgc_specimen_id','exp_value']
    df_probe = df_probe.drop_duplicates(subset=['icgc_specimen_id'])
    df_gene = df_gene.drop_duplicates(subset=['icgc_specimen_id'])
    df_concat = df_probe.join(df_gene.set_index('icgc_specimen_id'), on='icgc_specimen_id')
    low_meth = df_concat[df_concat.meth_value.argsort().argsort() < len(df_concat)*0.2]
    high_meth = df_concat[df_concat.meth_value.argsort().argsort() > len(df_concat)*0.8]
    stats, pvalue = mwu(low_meth.exp_value, high_meth.exp_value)
    return stats, pvalue

##### load prot
df_prot = meth_prot()
##### load ensembl gene
df_ensembl = ensembl_gene()

##### recurrently use MWU test to detect significant correlation between probe and gene
proj_code ='BRCA-US'
i_th = 0
df_meth = pd.read_hdf('./data/meth_matrix_hdf/{}/meth_matrix.{}_{}.hdf'.format(proj_code, proj_code, i_th))
df_exp = pd.read_hdf('./data/exp_matrix/exp_matrix.{}.hdf'.format(proj_code))
result = []
for probe_id in probe_list:
    for gene_id in gene_list:
        probe_chr, probe_pos = df_prot[['CHR','MAPINFO']][df_prot.IlmnID==probe_id].values[0]
        gene_chr, gene_start, gene_end, gene_strand = df_ensembl[['chr','start','end','strand']][df_ensembl.gene_id==gene_id].values[0]
        df_probe = df_meth[probe_id]
        df_gene = df_exp[gene_id]
        # WM U-test
        stats, pvalue = eQTM(df_probe, df_gene, test='mwu')
        result.append([probe_id, probe_chr, probe_pos, gene_id, gene_chr, gene_start, gene_end, gene_strand, stats, pvalue])
##### transfer all result into a dataframe
df_result = pd.DataFrame(result, columns=['probe','probe_chr','probe_pos','gene','gene_chr','gene_start','gene_end','gene_strand','u-stat','pvalue'])