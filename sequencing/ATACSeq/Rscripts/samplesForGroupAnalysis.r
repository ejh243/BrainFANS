## ===============================================================================================================##
##                                PIPPA STEP 8.0: Sample list for group peak calling                              ##
## ===============================================================================================================##
## EXECUTION: Rscript ./Rscripts/samplesForGroupAnalysis.r <project_config> <cf>                                  ||
##                                                                                                                ||
## DESCRIPTION: This scripts creates a file with samples that passed both stage 1 and 2 of the QC and selects     ||
## samples that belong to an input cell type to perform peak calling.                                             ||
##                                                                                                                ||
## INPUTS:                                                                                                        ||
## - <project_config> : path to config.r file.                                                                    ||
## - <cf>: cell fraction of samples to select for group peak calling                                              ||
##                                                                                                                ||
## OUTPUTS:                                                                                                       ||
## - "samplesForGroupAnalysisOrdered_<cf>.txt" in 0_metadata folder                                               ||
## - "samplesGroupAnalysis.csv" in 0_metadata folder                                                              ||
##                                                                                                                ||
## REQUIRES:                                                                                                      ||
## - R version > 4.3                                                                                             ||
## - A csv file in the metadata folder with the summary of QC stages 1 and 2: passS1S2Status.csv. This is produced||
##   at STEP 7.2 SUMMARY                                                                                          ||
## ===============================================================================================================##

args <- commandArgs(trailingOnly=TRUE)
configFile <-args[1]
source(configFile)
cf <- args[2]

## collate only samples that passed Stage 1 and 2 of QC
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
passAllQC <- read.csv(file.path(qc12stats), stringsAsFactors = FALSE, strip.white = TRUE)
passAllQC <- passAllQC[passAllQC$PASSALL == TRUE ,]
phenoAllQC <- pheno[pheno$sampleCode %in% passAllQC$sampleCode,]

print(table(phenoAllQC$fraction))
print(table(phenoAllQC$cohort))

## create csv file with only samples that passed Stage 1 and 2 of QC
write.csv(phenoAllQC, file = file.path(groupAnalysisSamples), row.names = FALSE)

passAllQC <- read.csv(file.path(groupAnalysisSamples), stringsAsFactors = FALSE, strip.white = TRUE)
samples <- passAllQC[passAllQC$fraction == cf,]$sampleID
samples <- na.omit(samples)

## output lits of samples for chosen cell type
write.table(samples, file = file.path(paste0(metaDir, "/samplesForGroupAnalysis_",cf,".txt")), 
            row.names = FALSE,quote = FALSE, col.names=FALSE)