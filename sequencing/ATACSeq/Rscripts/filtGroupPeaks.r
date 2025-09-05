## ===============================================================================================================##
##                          ATAC-seq pipeline STEP 7.2: Filter peaks called in groups                             ##
## ===============================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/filtGroupPeaks.r <project> <cell-group>                       ||
## - execute from scripts directory                                                                               ||
##                                                                                                                ||
## DESCRIPTION: This scripts filters the peaks called at group level based on whether peaks are found in at least ||
## one sample.                                                                                                    ||
##                                                                                                                ||
## INPUTS:                                                                                                        ||
## - <project> : project on which analysis is being run                                                           ||
## - <ell-group>: cell fraction of samples to select for group peak calling                                       ||
##                                                                                                                ||
## OUTPUTS:                                                                                                       ||
## - <cell-group>.filt.broadPeak.bed" in 4_calledpeaks folder                                                     ||
##                                                                                                                ||
## REQUIRES:                                                                                                      ||
## - R version > 4.3                                                                                              ||
## - A txt file with a list of cell-type samples samplesForGroupAnalysisOrdered_<cell-type>.txt in metadata folder//
## - Parameters in config.r file: samplesPeaks, sampleSheet, groupsPeaks,metaDir, macs_thres_group                ||
## ===============================================================================================================##



args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
cf<-args[2]
source(configFile)

library(GenomicRanges)
library(UpSetR)
library(rtracklayer)
library(dplyr)

ReadPeakstoGRanges <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand",
                    "signalValue", "pValue", "qValue")
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, end = df$end))
  return(gr)
}

## Gather samples peak sets
peak_files <- list.files(path = samplesPeaks, pattern = paste0("*\\.narrowPeak.filt$"), full.names = TRUE)

## Select samples for group analysis
samples <- read.table(file.path(paste0(metaDir, "/samplesForGroupAnalysisOrdered_",cf,".txt")))[,1]
pheno<-read.csv(file.path(sampleSheet))
pheno <- pheno[pheno$sampleID %in% samples,]

# Extract file names without extensions for labeling
file_labels <- tools::file_path_sans_ext(basename(peak_files))
file_labels <- unlist(lapply(file_labels, function(x) gsub(".narrowPeak.filt","",x)))
file_labels <- unlist(lapply(file_labels, function(x) gsub(".narrowPeak","",x)))
print("Number of samples that passed previous QC stages found: ")
nrow(pheno)
names(peak_files) <- file_labels
peak_files <- peak_files[names(peak_files) %in% pheno$sampleID]
print("Names of samples peaks found: ")
print(names(peak_files))
peak_sets <- lapply(peak_files, ReadPeakstoGRanges)

## Read cell-group peaks
peak_file_group <- paste0(groupsPeaks,cf,"_", macs_thres_group,".sorted.narrowPeak.filt")
peak_set_group <- ReadPeakstoGRanges(peak_file_group)
names(peak_set_group) <-  paste0(seqnames(peak_set_group), "-", start(peak_set_group), "-", end(peak_set_group))

## Find overlaps and create a binary matrix
overlap_matrix <- data.frame(matrix(0, nrow = length(peak_set_group), ncol = length(peak_sets)))
colnames(overlap_matrix) <- names(peak_sets)
for (i in seq_along(peak_sets)) {
  overlaps <- findOverlaps(peak_set_group, peak_sets[[i]])
  overlap_matrix[unique(queryHits(overlaps)), i] <- 1
}
rownames(overlap_matrix) <- names(peak_set_group)
print(head(names(peak_set_group)))
filtered_mat <- overlap_matrix[rowSums(overlap_matrix) > 0, ]
matrix_names <- rownames(filtered_mat)

# Filter GRanges based on matrix row names
filtered_gr <- peak_set_group[names(peak_set_group) %in% matrix_names]
export(filtered_gr, con = paste0(groupsPeaks,cf,"_", macs_thres_group,".filt.narrowPeak.bed"), format = "bed")
