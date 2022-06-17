##---------------------------------------------------------------------#
##
## Title: Calculate ChIP-seq QC metrics
##
## Purpose of script: to output sample summary QC metrics according to the CHIPQC r package
##
## Author: Jessica Shields
##
## Date Created: 2022-03-22
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## load arguments
args = commandArgs(trailingOnly=TRUE)
args[1]<-"epiGaba"
batchNum<-0

## load config variables
project<-args[1]
batchNum<-as.numeric(args[2]) ## nb starts from 0
source("ChIPSeq/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library(magrittr)
register(DoparParam())
registered() 
bpparam("SerialParam")

#----------------------------------------------------------------------#
# IMPORT AND WRANGLE DATA
#----------------------------------------------------------------------#
## Create sample sheet for the chipqc input
if (file.exists(paste0(metaDir, "/sampleSheetForChipQC.csv"))==FALSE){
  sampleSheet<-read.csv(sampleSheet)
  
  peaks<-list.files(peakDir, pattern = "Peak.filt", recursive = TRUE) %>%
    sort()
  bamReads<-paste(sampleSheet$sampleID, 'filt.nodup.bam') %>%
    sort()
  
  # bam files
  bamControl<-paste(sampleSheet$controlID, 'filt.nodup.bam') %>%
    sort()
  
  # necessary columns
  factor<-sampleSheet$target
  tissue<-sampleSheet$fraction
  pe<-"Paired"
  peakIndex<-match(sampleSheet$sampleID, gsub(".narrowPeak.filt|.broadPeak.filt", "", peaks))
  
  sampleSheet<-data.frame(SampleID = sampleSheet$sampleID, Tissue=tissue, Factor=factor, Replicate=1, ReadType = pe, 
                          bamReads = paste(alignedDir, bamReads, sep = "/"), 
                          ControlID = sampleSheet$controlID,
                          bamControl = paste(alignedDir, bamControl, sep = "/"),
                          Peaks = paste(peakDir, peaks[peakIndex],sep = "/"),
                          PeakCaller='macs',
                          stringsAsFactors = FALSE)
  
  if (batchNum == 0){
    write.csv(sampleSheet, paste(metaDir, "sampleSheetForChipQC.csv",sep = "/"), row.names = FALSE)
  } 
} else {
  print('Using existing sampleSheet for ChIPQC')
  sampleSheet<- read.csv(paste0(metaDir,"/sampleSheetForChipQC.csv"))
}

#----------------------------------------------------------------------#
# SUBSET TO BATCH SUBMIT
#----------------------------------------------------------------------#
## filter to subset of samples
index<-c(1:10)+(batchNum*10)
## if number of samples is not a function of ten adjust index
index<-index[index %in% 1:length(rownames(sampleSheet))]

sampleSheet<- sampleSheet[index,]

#----------------------------------------------------------------------#
# CREATE CHIP QC OBJECT
#----------------------------------------------------------------------#

dat<-ChIPQC(sampleSheet, consensus = TRUE, chromosomes = NULL, annotation = "hg38", blacklist = blacklist)
save(dat, file = paste0(peakDir, "/QCOutput/ChIPQCObject_", batchNum, ".rdata"))



#sampleSheet$Peaks<-paste0(peakDir, "/", sampleSheet$Tissue, "_peaks.narrowPeak")
#dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")
#save(dat, file = paste(dataDir, "ChIPQCObjectFractionPeaks.rdata", sep = "/"))
