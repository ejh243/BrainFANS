## Title: Run cis MatrixEQTL interaction model
##
## Purpose of script: perform cis QTL analysis by chr across all samples 
## with interaction between genotype and cell type
##
## Author: Eilis Hannon
##
## Date Created: 2022-08-12
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(MatrixEQTL)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)
dataDir<-args[1]
chr<-args[2]
covFile<-args[3]

useModel = modelLINEAR_CROSS; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 
pvOutputThreshold = 1;
cisDist<-1000000
## set error covarriance. NOTE rarely used instead set to numeric()
errorCovariance = numeric();
 
SNP_file_name = file.path(dataDir, "genotype", paste0("genotype_", chr, ".txt"))
snps_location_file_name = file.path(dataDir, "genotype", paste0("genotype_", chr, "_Info.txt"))

methylation_file_name<-file.path(dataDir, "methylation", paste0("methylation_", chr, ".txt"))
methylation_location_file_name<-file.path(dataDir, "methylation", paste0("methylation_", chr, "_Info.txt"))

covariates_file_name = file.path(dataDir, "covariate", covFile) 
intVar<-gsub("\\.txt", "", tail(unlist(strsplit(covFile, "_")), n = 1))

output_file_name = file.path(sub("Input", "Output", dataDir), "InteractionModel", paste0("cisQTLs_chr", chr, "_CTInteraction_", intVar,".txt"))

if (!file.exists(file.path(sub("Input", "Output", dataDir), "InteractionModel"))) {
 dir.create(file.path(sub("Input", "Output", dataDir), "InteractionModel"), recursive = TRUE)
}

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);

## load SNP data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

## load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( methylation_file_name);

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
methpos = read.table(methylation_location_file_name, header = TRUE, stringsAsFactors = FALSE);


#----------------------------------------------------------------------#
# Run mQTLs
#----------------------------------------------------------------------#


me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	pvOutputThreshold = 0,
	output_file_name.cis = output_file_name,
	pvOutputThreshold.cis = pvOutputThreshold,
	snpspos = snpspos, 
	genepos = methpos,
	cisDist = cisDist,
	useModel = useModel, 
	errorCovariance = errorCovariance, 
	verbose = TRUE,
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE)
