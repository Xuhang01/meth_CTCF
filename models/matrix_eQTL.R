# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.
# command: for i in {0..276}; do echo $i; done | parallel -j10  Rscript ./models/matrix_eQTL.R
args <- commandArgs(T)
# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = paste("./data/",args[1],sep="");
# base.dir = paste("./data", args[1])
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/meth_matrix_ith/meth_matrix_", args[2], ".txt", sep="");
# SNP_file_name = args[2]
snps_location_file_name = paste(base.dir, "/meth450_probe.txt", sep="");
# snps_location_file_name = args[3]

# Gene expression file name
expression_file_name = paste(base.dir, "/exp_matrix.txt", sep="");
# expression_file_name = args[4]
gene_location_file_name = paste(base.dir, "/ensembl_gene.txt", sep="");
# gene_location_file_name = args[5]


# Covariates file name
# Set to character() for no covariates
# covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
covariates_file_name = c()

# Output file name
output_file_name_cis = paste(base.dir, "/results/", args[2] , "_cis_eQTL.txt", sep="");
output_file_name_tra = paste(base.dir, "/results/", args[2] , "_trans_eQTL.txt", sep="");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e4;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

##### show results
# unlink(output_file_name_tra);
# unlink(output_file_name_cis);

# ## Results:

# cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
# cat('Detected local eQTLs:', '\n');
# show(me$cis$eqtls)
# cat('Detected distant eQTLs:', '\n');
# show(me$trans$eqtls)

# ## Plot the Q-Q plot of local and distant p-values

# plot(me)
