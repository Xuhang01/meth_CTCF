import pandas  as pd
import numpy as np
import os
##### constants
proj_list = ['BLCA-US','BRCA-US','CESC-US','COAD-US','GBM-US','KIRC-US','KIRP-US','LAML-US',
             'LIHC-US','LUAD-US','LUSC-US','OV-US','PAAD-US','PRAD-US','READ-US','SKCM-US',
             'STAD-US','THCA-US','UCEC-US']
##### load methylation file
def load_meth(proj_code, i_th):
    return pd.read_hdf('./data/meth_matrix_hdf/{}/meth_matrix.{}_{}.hdf'.format(proj_code, proj_code, i_th))
def load_exp(proj_code):
    return pd.read_hdf('./data/exp_matrix/exp_matrix.{}.hdf'.format(proj_code))
##### load icgc specimen file
def icgc_specimen(proj_code):
    df_sp = pd.read_table("/f/hang/PhD_Project/Data/Database/ICGC/ICGC_Project/{}/specimen.{}.tsv.gz".format(proj_code, proj_code))
    df_sp.loc[:,"type"] = df_sp.specimen_type.apply(
            lambda s: 
            "PTS" if s == "Primary tumour - solid tissue" 
            else "NB" if s == "Normal - blood derived"
            else "NT" if s == "Normal - tissue adjacent to primary"
            else "MT" if s == "Metastatic tumour - metastasis to distant location"
            else "RTS" if s == "Recurrent tumour - solid tissue"
            else "RTO" if s == "Recurrent tumour - other"
            else "PTBP" if s == "Primary tumour - blood derived (peripheral blood)"
            else "NO" if s == "Normal - other"
            else "ERROR")
    return df_sp
##### load icgc_donor_file
def icgc_donor(proj_code):
    return pd.read_table("/f/hang/PhD_Project/Data/Database/ICGC/ICGC_Project/{}/donor.{}.tsv.gz".format(proj_code, proj_code))
##### load methylation protocol files
def meth_prot():
    columns_type = {
        'IlmnID': 'object',
        'Name': 'object',
        'AddressA_ID': 'object',
        'AlleleA_ProbeSeq': 'object',
        'AddressB_ID': 'object',
        'AlleleB_ProbeSeq': 'object',
        'Infinium_Design_Type': 'category',
        'Next_Base': 'category',
        'Color_Channel': 'category',
        'Forward_Sequence': 'object',
        'Genome_Build': 'category',
        'CHR': 'category',
        'MAPINFO': 'object',
        'SourceSeq': 'object',
        'Chromosome_36': 'category',
        'Coordinate_36': 'object',
        'Strand': 'category',
        'Probe_SNPs': 'object',
        'Probe_SNPs_10': 'object',
        'Random_Loci': 'category',
        'Methyl27_Loci': 'category',
        'UCSC_RefGene_Name': 'object',
        'UCSC_RefGene_Accession': 'object',
        'UCSC_RefGene_Group': 'category',
        'UCSC_CpG_Islands_Name': 'object',
        'Relation_to_UCSC_CpG_Island': 'category',
        'Phantom': 'object',
        'DMR': 'category',
        'Enhancer': 'category',
        'HMM_Island': 'object',
        'Regulatory_Feature_Name': 'object',
        'Regulatory_Feature_Group': 'category',
        'DHS': 'category'
    }
    meth_prot = '/f/hang/PhD_Project/Data/Database/ICGC/Methylation_Protocol_file/HumanMethylation450_15017482_v1-2.csv'
    df_prot = pd.read_csv(meth_prot, skiprows=range(7), dtype=columns_type)
    df_prot = df_prot[df_prot.IlmnID.apply(lambda s: s.startswith('cg'))]
    df_prot.loc[:,"order"] = df_prot.IlmnID.apply(lambda s: int(s[2:])//100000)
    df_prot = df_prot.sort_values(by=["IlmnID"])
    df_prot.MAPINFO = df_prot.MAPINFO.astype(np.int)
    return df_prot[['IlmnID','CHR','MAPINFO','order']]
##### load ensembl gene file
def ensembl_gene():
    ensembl_gene = "./data/Homo_sapiens.GRCh37.87.gtf"
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
    return df_gene
##### given probe_id, find the nearest genes (10)
def probe_nearest_gene(df_ensembl, df_prot, probe_id):
    # detect the most nearest genes: 10 upstream and 10 downstream
    probe_chr, probe_pos = df_prot[['CHR','MAPINFO']][df_prot.IlmnID == probe_id].values[0]
    sub_ensembl = df_ensembl[df_ensembl.chr == str(probe_chr)]
    sub_ensembl.loc[:,'distance'] = sub_ensembl.start - probe_pos
    upstream = sub_ensembl[sub_ensembl.distance < 0]
    downstream = sub_ensembl[sub_ensembl.distance > 0]
    ten_up = upstream.sort_values('distance')[-10:]
    ten_down = downstream.sort_values('distance')[:10]
    ten_up.loc[:,'pos_type'] = 'upstream'
    ten_down.loc[:,'pos_type'] = 'downstream'
    return pd.concat([ten_up, ten_down])

##### find donors with both Normal sample and Tumor Sample in meth and exp dataset
def common_donor(df_meth, df_exp, proj_code):
    ##### load specimen file
    df_sp = icgc_specimen(proj_code)
    ##### load donor file
    df_donor = icgc_donor(proj_code)
    df_sp.loc[:,'donor_sex'] = df_sp.icgc_donor_id.apply(lambda s: df_donor[df_donor.icgc_donor_id==s].donor_sex.iat[0])
    df_sp.loc[:,'donor_age'] = df_sp.icgc_donor_id.apply(lambda s: df_donor[df_donor.icgc_donor_id==s].donor_age_at_diagnosis.iat[0])
    ##### find donor with both Normal Tissue and Primary Tumor Solid Tissue Specimen
    total_donor = df_sp[['project_code','icgc_donor_id','icgc_specimen_id','type']]
    # methylation with both NT and PTS 
    meth_donor = df_meth[['project_code','icgc_donor_id','icgc_specimen_id']]
    meth_donor = meth_donor[['icgc_specimen_id']].join(total_donor.set_index('icgc_specimen_id'),on='icgc_specimen_id')
    if 'NT' not in meth_donor.type.unique() or 'PTS' not in meth_donor.type.unique():
        overlap = pd.DataFrame(columns=['icgc_donor_id','NT','PTS','donor_sex','donor_age'])
        overlap.to_csv('./data/common_donor/{}.common.txt'.format(proj_code),sep='\t',index=False)
        return False
    meth_donor_count = meth_donor[['icgc_donor_id','icgc_specimen_id','type']].drop_duplicates()
    df_meth_donor = pd.pivot_table(meth_donor_count, columns='type',
        index="icgc_donor_id",values='icgc_specimen_id',aggfunc='first')[['NT','PTS']].dropna(subset=['NT','PTS'])
    ##### exp with both NT and PTS
    exp_donor = df_exp[['project_code','icgc_donor_id','icgc_specimen_id']]
    exp_donor = exp_donor[['icgc_specimen_id']].join(total_donor.set_index('icgc_specimen_id'),on='icgc_specimen_id')
    if 'NT' not in exp_donor.type.unique() or 'PTS' not in exp_donor.type.unique():
        overlap = pd.DataFrame(columns=['icgc_donor_id','NT','PTS','donor_sex','donor_age'])
        overlap.to_csv('./data/common_donor/{}.common.txt'.format(proj_code),sep='\t',index=False)
        return False
    exp_donor_count = exp_donor[['icgc_donor_id','icgc_specimen_id','type']].drop_duplicates()
    df_exp_donor = pd.pivot_table(exp_donor_count, columns='type',
        index="icgc_donor_id",values='icgc_specimen_id',aggfunc='first')[['NT','PTS']].dropna(subset=['NT','PTS'])
    ##### overlap between meth and exp
    df_meth_donor.columns = ['NT_meth','PTS_meth']
    df_exp_donor.columns = ['NT_exp','PTS_exp']
    overlap = df_exp_donor.reset_index().join(df_meth_donor,on='icgc_donor_id')
    overlap = overlap.dropna(subset=['NT_meth','PTS_meth','NT_exp','PTS_exp'])
    print('There are {} donor with NT,PTS,meth,exp'.format(overlap.shape[0]))
    overlap = overlap[['icgc_donor_id','NT_meth','PTS_meth']]
    overlap.columns = ['icgc_donor_id','NT','PTS']
    overlap.loc[:,'donor_sex'] = overlap.icgc_donor_id.apply(lambda s: df_donor[df_donor.icgc_donor_id==s].donor_sex.iat[0])
    overlap.loc[:,'donor_age'] = overlap.icgc_donor_id.apply(lambda s: df_donor[df_donor.icgc_donor_id==s].donor_age_at_diagnosis.iat[0]) 
    overlap.to_csv('./data/common_donor/{}.common.txt'.format(proj_code),sep='\t',index=False)


def main():
    ##### find donors with both Normal sample and Tumor Sample
    # This is not very good because such sample is very few
    df_meth = pd.read_hdf('./data/meth_matrix_hdf/ICGC/meth_matrix.ICGC_0.hdf')
    df_exp = pd.read_hdf('./data/exp_matrix/exp_matrix.ICGC.hdf')
    df_meth = df_meth.reset_index()
    df_exp = df_exp.reset_index()
    for proj_code in proj_list:
        print(proj_code)
        sub_meth= df_meth[df_meth.project_code == proj_code]
        sub_exp = df_exp[df_exp.project_code == proj_code]
        print(sub_meth.shape)
        print(sub_exp.shape)
        common_donor(sub_meth, sub_exp, proj_code)
    ##### 