## ================================================================================================##
##                          ATAC-seq pipeline STEP 7.3: Counts in peaks                            ##
## ================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/countsInPeaks.r <project>                      ||
## - execute from scripts directory                                                                ||
##                                                                                                 ||
## DESCRIPTION: This scripts gets read counts in peaks called at cell group level. Only samples    ||
## that passed both QC stages are used, as done to call peaks.                                     ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <project> : project on which analysis is being run                                            ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## -  counts in peaks (all peaks) file: peakCounts.rdata in 5_countPeaks folder                    ||
## -  counts in peaks (only promoters) file: peakCounts_prom.rdata in 5_countPeaks folder          ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R/4.2.1-foss-2022a, R libraries: Rsubread, GenomicRanges, ChIPpeakAnno,ChIPseeker             ||
##   TxDb.Hsapiens.UCSC.hg38.knownGene,org.Hs.eg.db                                                ||
## - A csv file in the metadata folder with samples used for group peak calling:                   || 
##   samplesGroupAnalysis.csv. This is produced in STEP 7.0                                        ||
## ================================================================================================##

## ==========##
##   SET-UP  ##
## ==========##

args <- commandArgs()
configFile<-args[6]
source(configFile)

library(Rsubread)
library(GenomicRanges)
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

## ===============##
##   GATHER DATA  ##
## ===============##

peakDir<- paste0(dir, "/4_calledPeaks/MACS/BAMPE/group/")
samples <- read.csv(file.path(paste0(metaDir,"/samplesGroupAnalysis.csv")), sep = ",",stringsAsFactors = FALSE)
pheno<-read.csv(file.path(metaDir, "sampleSheet.csv"))
pheno <- pheno[pheno$sampleID %in% samples$sampleID,]

peakFiles<-list.files(peakDir, pattern="*.sorted.broadPeak.filt")
print("Peaks files: ")
print(peakFiles)
peaks<-list()
for(each in peakFiles){
  fraction<-head(unlist(strsplit(each, "\\.")), n = 1)
  peaks[[fraction]]<-makeGRangesFromDataFrame(read.table(file.path(peakDir, each)), keep.extra.columns=TRUE, seqnames.field="V1",
                                              start.field="V2",
                                              end.field="V3",
                                              ignore.strand=TRUE)
}
names(peaks) <- cellTypes

## Annotate peaks
peaksAnno<-lapply(peaks,  annotatePeak, tssRegion=c(-1000, 1000),TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", verbose = FALSE)

bamFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam$", recursive = TRUE, include.dirs = TRUE)
sampleNames<-gsub(".filt.nodup.bam$", "", bamFiles)
sampleNames <- sampleNames[sampleNames %in% pheno$sampleID]
bamFiles<-bamFiles[bamFiles %in% unlist(lapply(sampleNames, function(x) paste0(x,".filt.nodup.bam")))]

## ======================##
##   COUNTS (PROMOTERS)  ##
## ======================##

ctPromotorPeaks <- list()
for(i in 1:length(cellTypes)){
  ct <- cellTypes[i]
  ctPromotorPeaks[[i]] <- peaksAnno[[ct]]@anno[grep("Promoter", peaksAnno[[ct]]@anno$annotation),]
}
ctPromotorPeaks <- do.call(c, ctPromotorPeaks)
ctPeaks<-data.frame("GeneID"=ctPromotorPeaks$V4, "Chr"=seqnames(ctPromotorPeaks), "Start"=start(ctPromotorPeaks), "End"=end(ctPromotorPeaks), "Strand"="*")

## Get counts in peaks annotated as promoters
fc_ctPeaks <- featureCounts(file.path(paste0(alignedDir,"/", bamFiles)),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = paste0(dir,"/5_countPeaks/Counts/peakCounts_prom.rdata"))


## ================##
##   COUNTS (ALL)  ##
## ================##

allPeaks <- list()
for(i in 1:length(cellTypes)){
  ct <- cellTypes[i]
  allPeaks[[i]] <- peaksAnno[[ct]]@anno
}
allPeaks <- do.call(c, allPeaks)
ctPeaks<-data.frame("GeneID"=allPeaks$V4, "Chr"=seqnames(allPeaks), "Start"=start(allPeaks), "End"=end(allPeaks), "Strand"="*")

## Get counts in all peaks
fc_ctPeaks <- featureCounts(file.path(paste0(alignedDir,"/", bamFiles)),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = paste0(dir,"/5_countPeaks/Counts/peakCounts_all.rdata"))