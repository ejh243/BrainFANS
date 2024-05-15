

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
req_cols <- c('Sample_ID','Individual_ID','Chip_ID','Chip_Location') 	# required column names
opt_cols <- c('Age')													# optional column names
cond_cols <- c('Sex','Genotype_IID','Cell_Type')						# conditionally column names

# source checkColnames()
source("checkColnamesFunction.r")


#1. Check required column names ------------------------------------------------------------------------------
print("1. Checking required column names:")
print(req_cols)
checkColnames(sampleSheet, req_cols, type='Required')


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