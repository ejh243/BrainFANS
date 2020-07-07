## this script is to be run after peak calling as assesses aligned reads

source("config.r") ## contains paths to files

library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel") 
library(ATACseqQC)
library(Rsamtools)
library("phastCons100way.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38")
library(MotifDb)
register(DoparParam())
registered() 
bpparam("SerialParam")

setwd(dataDir) ## change to directory where aligned files (bam) and peaks (narrowPeaks) can be found ## will search for all within this folder

### Create Sample Sheet
Peaks<-list.files(peakDir, pattern = ".broadPeak", recursive = TRUE)
peakIDs<-gsub("_depDup_q30.bam_peaks.broadPeak", "", Peaks)
bamReads<-list.files(alignedDir, pattern = "_depDup_q30.bam", recursive = TRUE)
bamReads<-bamReads[endsWith(bamReads, "_depDup_q30.bam")]
bamIDs<-gsub("_depDup_q30.bam", "", bamReads)


sampleIDs<-intersect(bamIDs, peakIDs)

tissue<-substr(unlist(lapply(strsplit(sampleIDs, "-"),tail, n = 1)), 1,1)

pe<-"Paired"

sampleSheet<-data.frame(SampleID = sampleIDs, Tissue=tissue, Factor="ATAC", Replicate=1, 
ReadType = pe, bamReads = paste(alignedDir, bamReads[match(sampleIDs, bamIDs)], sep = "/"), 
Peaks = paste(peakDir, Peaks[match(sampleIDs, peakIDs)], sep = "/"), stringsAsFactors = FALSE)

write.csv(sampleSheet, "SampleSheetForATACSeqQC.csv") ## overwrites existing samples sheet

dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")

save(dat, file = "ATACSeqQCObject_2.rdata")

