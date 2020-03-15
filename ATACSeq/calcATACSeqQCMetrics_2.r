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
Peaks<-list.files(".", pattern = ".narrowPeak", recursive = TRUE)
#Peaks<-Peaks[startsWith(Peaks, "H")]
bamReads<-list.files(".", pattern = "_depDup_q30.bam", recursive = TRUE)
#bamReads<-bamReads[startsWith(bamReads, "H")]
bamReads<-bamReads[grep("bai", bamReads, invert = TRUE)]
bamIDs<-unlist(lapply(strsplit(bamReads, "/"), tail, n = 1))
bamIDs<-gsub("_trimmed_depDup_q30.bam", "", bamIDs)
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


tissue<-unlist(lapply(strsplit(peakIDs, "_"),tail, n = 1))

pe<-"Paired"

sampleSheet<-data.frame(SampleID = indexS, Tissue=tissue[indexP], Factor="ATAC", Replicate=repN, ReadType = pe, bamReads = bamReads[indexB], Peaks = Peaks[indexP], stringsAsFactors = FALSE)

write.csv(sampleSheet, "SampleSheetForATACSeqQC.csv") ## overwrites existing samples sheet

dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")

save(dat, file = "ATACSeqQCObject_2.rdata")

