

## This script checks that sample sheet columns are formatted correctly prior to DNAm QC ##


args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
configFile <- paste0(dataDir, "/config.r")

# Load libraries
library(stringdist) # for amatch()
'%ni%' <- Negate('%in%') # define '%ni%' (not in)

# Load sample sheet
sampleSheet <- read.csv(paste0(dataDir, "/0_metadata/sampleSheet.csv"), na.strings = c("", "NA"), stringsAsFactors = FALSE)

# Column names to test
req_cols <- c('Sample_ID','Individual_ID') 												# required column names
bsnm_cols <- c('Basename','Chip_ID','Chip_Location','Sentrix_ID','Sentrix_Position')	# required when Basename not present
opt_cols <- c('Age')																	# optional column names
cond_cols <- c('Sex','Genotype_IID','Cell_Type')										# conditionally column names

# source checkColnames()
source("checkColnamesFunction.r")


#1. Check required column names ------------------------------------------------------------------------------
print("1. Checking required column names:")
print(c('Basename',req_cols))

# check Basename first
b <- checkColnames(sampleSheet, bsnm_cols[1], type='Required', verbose=F)

if(b$allPresent){
	# if Basename present, continue to check Basename alongside other required columns
	checkColnames(sampleSheet, c('Basename',req_cols), type='Required')
}else{
	# if Basename not present, check Chip and Sentrix as alternatives
	chip <- checkColnames(sampleSheet, bsnm_cols[2:3], type='Required', verbose=F)
	sntrx <- checkColnames(sampleSheet, bsnm_cols[4:5], type='Required', verbose=F)
	
	# if either Chip or Sentrix present, continue to check other required columns
	if(any(chip$allPresent | sntrx$allPresent)){
		cat("Basename column not found, but at least 2 of the following alternative columns are present: ", bsnm_cols[2:5])
		cat("Checking remaining required columns")
		checkColnames(sampleSheet, req_cols, type='Required')
	}else{
		cat(c("Basename column not found, and neither set of alternative column names are present: ", bsnm_cols[2:5]))
		cat("Checking remaining required columns")
		checkColnames(sampleSheet, req_cols, type='Required')
	}
}




#2. Check optional column names ------------------------------------------------------------------------------
print("2. Checking optional column names:")
print(opt_cols)
checkColnames(sampleSheet, opt_cols, type='Optional')


#3. Check conditional column names ---------------------------------------------------------------------------
print("Sourcing conditional variables from config.r")
source(configFile)
cond_status <- c(sexCheck,snpCheck,ctCheck) # T/Fs from config file
print(paste0(c("sexCheck","snpCheck","ctCheck"),"=", cond_status))

# subset conditional colnames to those TRUE in config
cond_cols.filtered <- cond_cols[cond_status]

if(all(cond_status==F)){
	print("No conditional variables to check")
}else{
	print("3. Checking conditional column names:")
	print(cond_cols.filtered)
	checkColnames(sampleSheet, cond_cols.filtered, type='Conditional')
}

# ---------------------------------------------------------------------------------------------------------- #