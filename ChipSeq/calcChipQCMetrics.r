
source("config.r") ## contains paths to files

library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel") 
register(DoparParam())
registered() 
bpparam("SerialParam")

setwd(dataDir) ## change to directory where aligned files (bam) and peaks (narrowPeaks) can be found ## will search for all within this folder

### Create Sample Sheet
Peaks<-list.files(peakDir, pattern = ".narrowPeak", recursive = TRUE)
bamReads<-list.files(alignedDir, pattern = "_sorted.bam", recursive = TRUE)
bamReads<-bamReads[grep("bai", bamReads, invert = TRUE)]
sampleIDs<-intersect(gsub("_peaks.narrowPeak", "", Peaks), gsub("_sorted.bam", "", bamReads))
tissue<-unlist(lapply(strsplit(sampleIDs, "_"), tail, n = 1))
pe<-"Paired"

peakIndex<-match(sampleIDs, gsub("_peaks.narrowPeak", "", Peaks))

sampleSheet<-data.frame(SampleID = sampleIDs, Tissue=tissue, Factor="H3K27ac", Replicate=1, ReadType = pe, bamReads = paste(alignedDir, bamReads, sep = "/"), Peaks = paste(peakDir, Peaks[peakIndex],sep = "/"), stringsAsFactors = FALSE)

write.csv(sampleSheet, paste(dataDir, "SampleSheetForChipQC.csv",sep = "/"))

dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")

save(dat, file = paste(dataDir, "ChIPQCObject.rdata", sep = "/"))

sampleSheet$Peaks<-paste0(peakDir, "/", sampleSheet$Tissue, "_peaks.narrowPeak")


dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")

save(dat, file = paste(dataDir, "ChIPQCObjectFractionPeaks.rdata", sep = "/"))
