## ================================================================================================##
##                             PIPPA STEP 10.3: Counts in curated peaks                            ##
## ================================================================================================##
## EXECUTION: Rscript ./Rscripts/countsInCuratedPeaks.r <project_config> <cf> <thres>              ||
##                                                                                                 ||
## DESCRIPTION: This scripts gets read counts in curated peaks called at chosen cell type. Only    ||
## samples that passed all QC stages are used, as done to call peaks. Peak set used will be all,   ||
## filtered for blacklisted regions only.                                                          ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <project_config> : path to config.r file.                                                     ||
## - <cf> : chosen cell type                                                                       ||
## - <thres> : q value threshold used for curated peak set                                         ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## -  counts in peaks (all peaks) file: in 5_countPeaks folder "downstreamCounts"                  ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R version > 4.3, R libraries: Rsubread, GenomicRanges                                         ||
## - A txt file with list of samples used for curated peak calling: downAnalysisSamples_<cf>.txt   ||
## ================================================================================================##

## ==========##
##   SET-UP  ##
## ==========##

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
cf<-args[2]
thres<-args[3]

source(configFile)

suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(GenomicRanges)
}))

## ===============##
##   GATHER DATA  ##
## ===============##

samples <- read.table(file.path(paste0(metaDir, "/downAnalysisSamples_",cf,".txt")))[,1]
pheno<-read.csv(file.path(sampleSheet))
pheno <- pheno[pheno$sampleID %in% samples,]

print(dim(pheno))
print(table(pheno$cohort))
print(table(pheno$fraction))

peakFile <- paste0(peaksCurated,"/", cf, "/", cf,"_",thres,".narrowPeak.filt")
print("Peaks file: ")
print(peakFile)
peaks<-makeGRangesFromDataFrame(read.table(peakFile), keep.extra.columns=TRUE, seqnames.field="V1",
                                              start.field="V2",
                                              end.field="V3",
                                              ignore.strand=TRUE)
                     
bamFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam$", recursive = TRUE, include.dirs = TRUE)
sampleNames<-gsub(".filt.nodup.bam$", "", bamFiles)
sampleNames <- sampleNames[sampleNames %in% samples]

bamFiles <- ifelse(
  grepl("_merged$", samples),
  paste0(samples, ".sorted.bam"),
  paste0(samples, ".filt.nodup.bam")
)

## ================##
##   COUNTS (ALL)  ##
## ================##

ctPeaks <- data.frame(peaks)
ctPeaks$GeneID <- paste0("Peak_", row.names(ctPeaks))
ctPeaks <- ctPeaks[c( "GeneID","seqnames","start","end","strand")]
names(ctPeaks) <- c("GeneID", "Chr", "Start", "End", "Strand")

## Get counts in all peaks
fc_ctPeaks <- featureCounts(file.path(paste0(alignedDir,"/", bamFiles)),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = paste0(downstreamCounts, "/peakCounts_",cf,"_", thres, ".narrowPeak.rdata"))