## ===============================================================================================================##
##                          ATAC-seq pipeline STEP 7.0: Sample list for group peak calling                        ##
## ===============================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/samplesForGroupAnalysis.r <project> <cell-group>              ||
## - execute from scripts directory                                                                               ||
##                                                                                                                ||
## DESCRIPTION: This scripts creates a file with samples that passed both stage 1 and 2 of the QC and selects     ||
## samples that belong to an input cell type to perform peak calling.                                             ||
##                                                                                                                ||
## INPUTS:                                                                                                        ||
## - <project> : project on which analysis is being run                                                           ||
## - <ell-group>: cell fraction of samples to select for group peak calling                                       ||
##                                                                                                                ||
## OUTPUTS:                                                                                                       ||
## - samplesForGroupAnalysisOrdered_<cell-group>.txt" in 0_metadata folder                                        ||
##                                                                                                                ||
## REQUIRES:                                                                                                      ||
## - R version > 4.3                                                                                              ||
## - A csv file in the metadata folder with the summary of QC stages 1 and 2: passS1S2Status.csv. This is produced||
##   at STEP 6.2 CHECK                                                                                            ||
## ===============================================================================================================##

args <- commandArgs(trailingOnly=TRUE)
configFile <-args[1]
source(configFile)
cf <- args[2]

pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
passAllQC <- read.csv(file.path(metaDir, "/passS1S2Status.csv"), stringsAsFactors = FALSE, strip.white = TRUE)
passAllQC <- passAllQC[match(pheno$sampleID,passAllQC$sampleID),]
passAllQC <- passAllQC[passAllQC$PASSALL == TRUE ,]
loSamples <- read.table(file.path(paste0(metaDir,"/leaveOutSamples.txt")))[,1]
passAllQC <- passAllQC[!passAllQC$sampleID %in% loSamples,]
write.csv(passAllQC, file = file.path(metaDir, "/samplesGroupAnalysis.csv"), row.names = FALSE)

samples <- passAllQC[passAllQC$fraction == cf,]$sampleID
samples <- na.omit(samples)
write.table(samples, file = file.path(paste0(metaDir, "/samplesForGroupAnalysisOrdered_",cf,".txt")), 
            row.names = FALSE,quote = FALSE, col.names=FALSE)
