## ================================================================================================##
##                              PIPPA STEP 10.2: Counts in curated cell type peaks                 ##
## ================================================================================================##
## EXECUTION: Rscript ./Rscripts/countsInCleanPeaksCT.r <project_config>                           ||
##                                                                                                 ||
## DESCRIPTION: This scripts gets read counts in curated peaks called from all cell types. Only    ||
## samples that passed all QC stages are used, as done to call peaks.                              ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <project_config> : path to config.r file.                                                     ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## - counts in peaks (all peaks) file: in 5_countPeaks folder "groupsCounts"                       ||
## - consensus peak set file: "consensusPeaks_5e-2_passAll.narrowPeak.bed" in 4_calledPeaks folder ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R > 4.3, R libraries: Rsubread, GenomicRanges, rtracklayer                                    ||
## - A csv file with samples used for curated peak calling: downAnalysisSamples.csv. This is       ||
## produced in STEP 9.0                                                                            ||
## - Parameters in config.r file: downAnalysisSamples, sampleSheet, peaksAnalysis, groupsCounts    ||
## macs_thres_group, alignedDir,groupCountsClean                                                   ||
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
  library(rtracklayer)
}))

## ===============##
##   GATHER DATA  ##
## ===============##

pheno <- read.csv(file.path(sampleSheet), sep = ",",stringsAsFactors = FALSE)
passSamples <- read.csv(file.path(downAnalysisSamples), sep = ",",stringsAsFactors = FALSE)
pheno <- pheno[pheno$sampleID %in% passSamples$sampleID,]
print("Number of samples that passed all QC stages are: ")
print(nrow(pheno))
print("Number of samples per cell-fraction")
print(table(pheno$fraction))

peakFiles <- list.files(
  path = peaksAnalysis,
  pattern = "narrowPeak\\.filt$",
  recursive = TRUE,
  full.names = TRUE
)

print(peakFiles)
peakFiles <- peakFiles[!grepl("\\.sorted\\.noSC\\.narrowPeak\\.filt$", peakFiles)]

print("Peaks files found are: ")
print(peakFiles)

peaks<-list()
for(each in peakFiles){
  fraction<-head(unlist(strsplit(each, "\\.")), n = 1)
  peaks[[fraction]]<-makeGRangesFromDataFrame(read.table(each), keep.extra.columns=TRUE, seqnames.field="V1",
                                              start.field="V2",
                                              end.field="V3",
                                              ignore.strand=TRUE)
}

peaks<-reduce(unlist(GRangesList(peaks)),min.gapwidth = 150)
export(peaks, con = paste0(peaksAnalysis,"/consensusPeaks_",macs_thres_group, "_passAll.narrowPeak.bed"), format = "bed")

bamFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam$", recursive = TRUE, include.dirs = TRUE)

sampleNames<-gsub(".filt.nodup.bam$", "", bamFiles)
sampleNames <- sampleNames[sampleNames %in% pheno$sampleID]

bamFiles<-bamFiles[bamFiles %in% unlist(lapply(sampleNames, function(x) paste0(x,".filt.nodup.bam")))]

## ================##
##   COUNTS (ALL)  ##
## ================##

## File name specified in config file
print("Counts in peaks file name is: ")
print(groupCountsClean)
           
allPeaks <- as.data.frame(peaks)
ctPeaks<-data.frame("GeneID"=paste0("Peak_", c(1:nrow(allPeaks))), "Chr"=allPeaks$seqnames, "Start"=allPeaks$start, "End"=allPeaks$end, "Strand"="*")
## Get counts in all peaks
fc_ctPeaks <- featureCounts(file.path(paste0(alignedDir,"/", bamFiles)),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = groupCountsClean)
print("Counts in peaks matrix done")

## ================##
##   DAR analysis  ##
## ================##
                                               
## Differential accessible analysis is done to get topmost differential peaks ##
## File name specified in config file
print("Differential peaks file name is: ")
print(diffPeaksCountsAll)
counts <- fc_ctPeaks$counts
colnames(counts) <- sub("\\.filt\\.nodup\\.bam$", "", colnames(counts))
cnt_table <- fc_ctPeaks$annotation
pheno.meta <- pheno[pheno$sampleID %in% colnames(counts),]
counts <- counts[,match(pheno.meta$sampleID,colnames(counts))]

dge <- DGEList(counts=counts,
           group=pheno.meta$fraction)
keep <- filterByExpr(dge$counts, group = pheno.meta$fraction)
dge <- dge[keep, , keep.lib.sizes = FALSE]

pheno.meta$fraction <- as.factor(pheno.meta$fraction)
print("Performing differential analysis to find top peaks")
## refCellType: base level cell-fraction should be specified in config file
pheno.meta$fraction = relevel(pheno.meta$fraction, ref = refCellType)
design <- model.matrix(~ fraction , data =pheno.meta )
d <- DGEList(counts=dge$counts,
           group=pheno.meta$fraction)
d <- calcNormFactors(d)
d <- estimateDisp(d, design, robust = TRUE)
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef=2) 
results <- lrt$table
results <- as.data.frame(results[order(results$PValue),])
results$Geneid <- rownames(results)
results$Bon <- p.adjust(results$PValue, method="bonferroni")
results$FDR <- p.adjust(results$PValue, method="fdr")
save(results,  file = diffPeaksCountsAll)