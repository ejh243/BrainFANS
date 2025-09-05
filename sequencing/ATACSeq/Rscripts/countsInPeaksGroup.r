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
## - <configFile> : config file of project on which analysis is being run                          ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## -  counts in peaks (all peaks) file: peakCounts.rdata in 5_countPeaks folder                    ||
## -  counts in peaks (only promoters) file: peakCounts_prom.rdata in 5_countPeaks folder          ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R version > 4.3, R libraries: Rsubread, GenomicRanges, ChIPpeakAnno,ChIPseeker                ||
##   TxDb.Hsapiens.UCSC.hg38.knownGene,org.Hs.eg.db                                                ||
## - A csv file in the metadata folder with samples used for group peak calling:                   || 
##   samplesGroupAnalysis.csv. This is produced in STEP 7.0                                        ||
## - Parameters in config.r file: groupAnalysisSamples, sampleSheet, macs_thres_group, groupsPeaks //
##   groupsCounts                                                                                  ||
## ================================================================================================##

## ==========##
##   SET-UP  ##
## ==========##

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(GenomicRanges)
  library(ChIPpeakAnno)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
}))

## ===============##
##   GATHER DATA  ##
## ===============##

samples <- read.csv(file.path(groupAnalysisSamples), sep = ",",stringsAsFactors = FALSE)
pheno<-read.csv(file.path(sampleSheet))
print(table(samples$fraction))
pheno <- pheno[pheno$sampleCode %in% samples$sampleCode,]
print(table(pheno$fraction))
peakFiles<-list.files(groupsPeaks, pattern=paste0("*_",macs_thres_group,".filt.narrowPeak.bed"))
print("Peaks files: ")
print(peakFiles)

peaks<-list()
for(each in peakFiles){
  fraction<-head(unlist(strsplit(each, "\\.")), n = 1)
  peaks[[fraction]]<-makeGRangesFromDataFrame(read.table(file.path(groupsPeaks, each)), keep.extra.columns=TRUE, seqnames.field="V1",
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
save(fc_ctPeaks, file = paste0(groupsCounts, "peakCounts_prom_filtbam_",macs_thres_group,".filt.narrowPeak.rdata"))


## ================##
##   COUNTS (ALL)  ##
## ================##

allPeaks <- list()
for(i in 1:length(names(peaks))){
  ct <- names(peaks)[i]
  allPeaks[[i]] <- peaksAnno[[ct]]@anno
}
allPeaks <- do.call(c, allPeaks)
ctPeaks<-data.frame("GeneID"=allPeaks$V4, "Chr"=seqnames(allPeaks), "Start"=start(allPeaks), "End"=end(allPeaks), "Strand"="*")

## Get counts in all peaks
fc_ctPeaks <- featureCounts(file.path(paste0(alignedDir,"/", bamFiles)),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = paste0(groupsCounts, "peakCounts_all_filtbam_",macs_thres_group,".filt.narrowPeak.rdata"))

