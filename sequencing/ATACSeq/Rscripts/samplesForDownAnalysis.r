## ===============================================================================================================##
##                         PIPPA STEP 9.0: Sample list for downstream analyses                                    ##
## ===============================================================================================================##
## EXECUTION: Rscript ./Rscripts/samplesForDownAnalysis.r <project_config> <cell-group>                           ||
##                                                                                                                ||
## DESCRIPTION: This scripts creates a file with samples that passed all stages of the QC and selects             ||
## samples that belong to an input cell type to perform peak calling.                                             ||
##                                                                                                                ||
## INPUTS:                                                                                                        ||
## - <project_config> : path to config.r file.                                                                    ||
## - <cell-group>: cell fraction of samples to select for group peak calling                                      ||
##                                                                                                                ||
## OUTPUTS:                                                                                                       ||
## - "downAnalysisSamples_<cell-group>.txt" in 0_metadata folder                                                  ||
## - "samplesDownAnalysis.csv" in in 0_metadata folder                                                            || 
##                                                                                                                ||
## REQUIRES:                                                                                                      ||
## - R > 4.3                                                                                                      ||
## - A csv file in the metadata folder with the summary of all QC stages: passAllQCStatus.csv This is             ||
## produced at STEP 8.4 CTCHECK                                                                                   ||
## ===============================================================================================================##

args <- commandArgs(trailingOnly=TRUE)
configFile <-args[1]
cf <- args[2]
source(configFile)

## collate only samples that passed all QC stages
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
passAllQC <- read.csv(file.path(allQCresults), stringsAsFactors = FALSE, strip.white = TRUE)
passAllQC <- passAllQC[passAllQC$PASSALL == TRUE ,]
phenoAllQC <- pheno[pheno$sampleCode %in% passAllQC$sampleCode,]

print(table(phenoAllQC$fraction))
print(table(phenoAllQC$cohort))

## create csv file with only samples that passed all QC
write.csv(phenoAllQC, file = file.path(downAnalysisSamples), row.names = FALSE)

samples <- phenoAllQC[phenoAllQC$fraction == cf, ]$sampleID
samples <- na.omit(samples)

## output lits of samples for chosen cell type
write.table(samples, file = file.path(paste0(metaDir, "/downAnalysisSamples_",cf,".txt")), row.names = FALSE,quote = FALSE, col.names=FALSE)