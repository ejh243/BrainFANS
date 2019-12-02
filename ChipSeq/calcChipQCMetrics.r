
source("config.r") ## contains paths to files

library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel") 
register(DoparParam())
registered() 
bpparam("SerialParam")

setwd(dataDir) ## change to directory where aligned files (bam) and peaks (narrowPeaks) can be found ## will search for all within this folder

### Create Sample Sheet
Peaks<-list.files(".", pattern = ".narrowPeak", recursive = TRUE)
bamReads<-list.files(".", pattern = "_sorted.bam", recursive = TRUE)
bamReads<-bamReads[grep("bai", bamReads, invert = TRUE)]
bamIDs<-unlist(lapply(strsplit(bamReads, "/"), tail, n = 1))
bamIDs<-gsub("_trimmed_sorted.bam", "", bamIDs)
peakIDs<-unlist(lapply(strsplit(Peaks, "/"), tail, n = 1))
peakIDs<-gsub("_trimmed_depDup_q30_peaks.narrowPeak", "", peakIDs)


sampleIDs<-intersect(bamIDs, peakIDs)
indexB<-NULL
indexP<-NULL
indexS<-NULL
repN<-NULL

for(each in sampleIDs){
	tmp.b<-which(bamIDs == each)
	tmp.p<-which(peakIDs == each)
	if(length(tmp.b) == length(tmp.p)){
		indexB<-c(indexB, tmp.b)
		indexP<-c(indexP, tmp.p)
		repN<-c(repN, 1:length(tmp.b))
		indexS<-c(indexS, rep(each, length(tmp.b)))
	} else {
		print(paste("Can't find matching BAM and PEAK files for", each))
	}
	
}


tissue<-"TOTAL"

pe<-"Paired"

sampleSheet<-data.frame(SampleID = indexS, Tissue=tissue, Factor="H3K27ac", Replicate=repN, ReadType = pe, bamReads = bamReads[indexB], Peaks = Peaks[indexP], stringsAsFactors = FALSE)

write.csv(sampleSheet, "SampleSheetForChipQC.csv")

dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")

save(dat, file = "ChIPQCObject.rdata")


