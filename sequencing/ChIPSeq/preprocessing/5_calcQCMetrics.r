##---------------------------------------------------------------------#
##
## Title: Calculate ChIP-seq QC metrics
##
## Purpose of script: to output sample summary QC metrics according to the CHIPQC r package
##
## Author: Eilis Hannon
##
## Date Created: 2022-03-22
##
##---------------------------------------------------------------------#

## load arguments
args = commandArgs(trailingOnly=TRUE)
<<<<<<< HEAD
=======
#args[1]<-"epiGaba"
#batchNum<-1
>>>>>>> e6d8709cbe265e69dd697458de440dd3dd2f6d2e

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
library(stringr)
register(DoparParam())
registered() 
bpparam("SerialParam")

#----------------------------------------------------------------------#
# IMPORT AND WRANGLE DATA
#----------------------------------------------------------------------#
## Create sample sheet
if (file.exists(paste0(metaDir, "/sampleSheetForChipQC.csv"))==FALSE){
  peaks<-list.files(peakDir, pattern = "Peak.filt", recursive = TRUE) %>%
    sort()
  bamReads<-list.files(alignedDir, pattern = "filt.nodup.bam", recursive = TRUE) %>%
    sort()
  bamReads<-bamReads[grep("bai", bamReads, invert = TRUE)]
  
  # bam files
  bamControl<-bamReads[grep("input", bamReads)]
  bamReads<- bamReads[grep("input", bamReads, invert=TRUE)]
  
  # create sample and control IDs 
  sampleIDs<-intersect(gsub(".narrowPeak.filt|.broadPeak.filt", "", peaks), gsub(".filt.nodup.bam", "", bamReads))
  controlIDs<- gsub(".filt.nodup.bam", "", bamControl)
  
  # necessary columns
  factor<-unlist(lapply(strsplit(sampleIDs, "_"), tail, n = 1)) %>%
    str_extract(., '\\b\\w+$') 
  tissue<-str_extract(sampleIDs, '\\.[A-Z]+') %>%
    sub('.', '', .) %>%
    str_replace(., 'SOX', 'GABA')
  pe<-"Paired"
  peakIndex<-match(sampleIDs, gsub(".narrowPeak.filt|.broadPeak.filt", "", peaks))
  
  sampleSheet<-data.frame(SampleID = sampleIDs, Tissue=tissue, Factor=factor, Replicate=1, ReadType = pe, 
                          bamReads = paste(alignedDir, bamReads, sep = "/"), 
                          ControlID = controlIDs,
                          bamControl = paste(alignedDir, bamControl, sep = "/"),
                          Peaks = paste(peakDir, peaks[peakIndex],sep = "/"),
                          PeakCaller='macs',
                          stringsAsFactors = FALSE)
  
  if (batchNum == 0){
    write.csv(sampleSheet, paste(metaDir, "sampleSheetForChipQC.csv",sep = "/"), row.names = FALSE)
  } 
} else if (file.exists(paste0(metaDir, "/sampleSheetForChipQC.csv"))==TRUE){
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
