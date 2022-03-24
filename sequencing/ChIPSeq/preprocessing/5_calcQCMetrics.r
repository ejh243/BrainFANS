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
#args[1]<-"epiGaba"

## load config variables
project<-args[1]
source("ChIPSeq/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel") 
library(magrittr)
library(stringr)
register(DoparParam())
registered() 
bpparam("SerialParam")

#----------------------------------------------------------------------#
# IMPORT AND WRANGLE DATA
#----------------------------------------------------------------------#

## Create sample sheet
peaks<-list.files(peakDir, pattern = ".narrowPeak.filt", recursive = TRUE)
bamReads<-list.files(alignedDir, pattern = "_sorted.bam", recursive = TRUE)
bamReads<-bamReads[grep("bai", bamReads, invert = TRUE)]

# control file and IDs
bamControl<-bamReads[grep("input", bamReads)]
controlIDs<- gsub("_depDup_q30.bam", "", bamControl)

# sample files and IDs
bamReads<- bamReads[grep("input", bamReads, invert=TRUE)]
sampleIDs<-intersect(gsub(".narrowPeak.filt", "", peaks), gsub("_sorted.bam", "", bamReads))

# necessary columns
factor<-unlist(lapply(strsplit(sampleIDs, "_"), tail, n = 1)) %>%
  str_extract(., '\\b\\w+$') 
tissue<-str_extract(sampleIDs, '\\.[A-Z]+') %>%
  sub('.', '', .)
pe<-"Paired"
peakIndex<-match(sampleIDs, gsub(".narrowPeak.filt", "", peaks))

sampleSheet<-data.frame(SampleID = sampleIDs, Tissue=tissue, Factor=factor, Replicate=1, ReadType = pe, 
                        bamReads = paste(alignedDir, bamReads, sep = "/"), 
                        ControlID = controlIDs,
                        bamControl = paste(alignedDir, bamControl, sep = "/"),
                        Peaks = paste(peakDir, peaks[peakIndex],sep = "/"), 
                        PeakCaller = 'macs',
                        stringsAsFactors = FALSE)

write.csv(sampleSheet, paste(metaDir, "sampleSheetForChipQC.csv",sep = "/"))


#----------------------------------------------------------------------#
# CREATE CHIP QC OBJECT
#----------------------------------------------------------------------#

dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")
save(dat, file = paste(peakDir, "QCOutput/ChIPQCObject.rdata", sep = "/"))



#sampleSheet$Peaks<-paste0(peakDir, "/", sampleSheet$Tissue, "_peaks.narrowPeak")
#dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")
#save(dat, file = paste(dataDir, "ChIPQCObjectFractionPeaks.rdata", sep = "/"))
