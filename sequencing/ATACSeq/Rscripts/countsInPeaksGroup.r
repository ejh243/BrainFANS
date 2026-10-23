## ================================================================================================##
##                               PIPPA STEP 8.3: Counts in peaks                                   ##
## ================================================================================================##
## EXECUTION: Rscript ./Rscripts/countsInPeaksGroup.r <project_config>                             ||
##                                                                                                 ||
## DESCRIPTION: This scripts gets read counts in peaks called at cell group level. Only samples    ||
## that passed both QC stages are used, as done to call peaks.                                     ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <project_config> : path to config.r file.                                                     ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## - consensus peak set across cell types: "consensusPeaks_5e-2_passS1-2.narrowPeak.bed" in        ||
## 4_calledPeaks folder                                                                            ||
## - counts in peaks (all peaks) file: peakCounts.rdata in 5_countPeaks folder                     ||
## - differential accessibility results file: diffPeaksCountsAll.rdata in 5_countPeaks folder      ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R version > 4.3, R libraries: Rsubread, GenomicRanges, edgeR, rtracklayer                     ||
## - A csv file in the metadata folder with samples used for group peak calling:                   || 
##   samplesGroupAnalysis.csv. This is produced in STEP 8.0                                        ||
## - Parameters in config.r file: groupAnalysisSamples, sampleSheet, macs_thres_group, groupsPeaks ||
##   groupCountsAll, diffPeaksCountsAll, alignedDir                                                ||
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
  library(edgeR)
}))

## ===============##
##   GATHER DATA  ##
## ===============##

pheno <- read.csv(file.path(sampleSheet), sep = ",",stringsAsFactors = FALSE)
passSamples <- read.csv(file.path(groupAnalysisSamples), sep = ",",stringsAsFactors = FALSE)
pheno <- pheno[pheno$sampleID %in% passSamples$sampleID,]
print("Number of samples that passed both stage 1 and 2 of QC is: ")
print(nrow(pheno))
print("Number of samples per cell-fraction")
print(table(pheno$fraction))

peakFiles <- list.files(
  path = groupsPeaks,
  pattern = "narrowPeak\\.filt$",
 recursive = TRUE,
  full.names = TRUE
)

peakFiles <- peakFiles[!grepl("\\.sorted\\.noSC\\.narrowPeak\\.filt$", peakFiles)]

print("Peaks files found are: ")
print(peakFiles)

peaks<-list()
for(each in peakFiles){
  fraction<-head(unlist(strsplit(each, "\\.")), n = 1)
  peaks[[fraction]]<-makeGRangesFromDataFrame(read.table(file.path(groupsPeaks, each)), keep.extra.columns=TRUE, seqnames.field="V1",
                                              start.field="V2",
                                              end.field="V3",
                                              ignore.strand=TRUE)
}

peaks<-reduce(unlist(GRangesList(peaks)),min.gapwidth = 150)
## save consensus peak set ##
export(peaks, con = paste0(groupsPeaks,"/consensusPeaks_",macs_thres_group,"_passS1-2.narrowPeak.bed"), format = "bed")

## gather bam files for read counts in peaks ##
bamFiles<-list.files(alignedDir, pattern = "*.filt.nodup.bam$", recursive = TRUE, include.dirs = TRUE)
sampleNames<-gsub(".filt.nodup.bam$", "", bamFiles)
sampleNames <- sampleNames[sampleNames %in% pheno$sampleID]

bamFiles<-bamFiles[bamFiles %in% unlist(lapply(sampleNames, function(x) paste0(x,".filt.nodup.bam")))]
stopifnot(all(dim(bamFiles) == nrow(pheno)))
bamFiles <- paste0(alignedDir,"/", bamFiles)

## ================##
##   COUNTS (ALL)  ##
## ================##

## File name specified in config file
print("Counts in peaks file name is: ")
print(groupCountsAll)

ctPeaks <- data.frame(peaks)
ctPeaks$GeneID <- paste0("Peak_", row.names(ctPeaks))
ctPeaks <- ctPeaks[c( "GeneID","seqnames","start","end","strand")]
names(ctPeaks) <- c("GeneID", "Chr", "Start", "End", "Strand")
print("Getting counts in peaks")
## Get counts in all peaks
fc_ctPeaks <- featureCounts(file.path(bamFiles),annot.ext=ctPeaks, allowMultiOverlap = TRUE, isPairedEnd = TRUE, nthreads = 10, fracOverlap=0.2)
save(fc_ctPeaks, file = groupCountsAll)
print("Counts in peaks matrix done")

## ================##
##   DAR analysis  ##
## ================##
                                               
## Differential accessible analysis is done to get topmost differential peaks ##
## File name specified in config file
print("Differential peaks file name is: ")
print(diffPeaksCountsClean)
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
save(results,  file = diffPeaksCountsClean)