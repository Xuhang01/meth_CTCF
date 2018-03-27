##### load methylation protocol files
def meth_prot():
    meth_prot = '/f/hang/PhD_Project/Data/Database/ICGC/Methylation_Protocol_file/HumanMethylation450_15017482_v1-2.csv'
    df_prot = pd.read_csv(meth_prot, skiprows=range(7))
    df_prot = df_prot[df_prot.IlmnID.apply(lambda s: s.startswith('cg'))]
    df_prot.loc[:,"order"] = df_prot.IlmnID.apply(lambda s: int(s[2:])//100000)
    df_prot = df_prot.sort_values(by=["IlmnID"])
    df_prot.MAPINFO = df_prot.MAPINFO.astype(np.int)
    return df_prot

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

def probe_nearest_gene(df_ensembl, df_prot, probe_id):
    # detect the most nearest genes: 10 upstream and 10 downstream
    probe_chr, probe_pos = df_prot[['CHR','MAPINFO']][df_prot.IlmnID == probe_id].values[0]
    sub_ensembl = df_ensembl[df_ensembl.chr == str(probe_chr)]
    sub_ensembl.loc[:,'distance'] = sub_ensembl.start - probe_pos

    upstream = sub_ensembl[sub_ensembl.distance < 0]
    downstream = sub_ensembl[sub_ensembl.distance > 0]

    ten_up = upstream.sort_values('distance')[-10:]
    ten_down = downstream.sort_values('distance')[:10]


