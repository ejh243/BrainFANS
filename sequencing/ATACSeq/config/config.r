## ===========================================================##
##           ATAC-seq pipeline: Example config R file         ##
## ===========================================================##
##                                                            ||
## This is an example config R file for the ATAC-sequencing   ||
## pipeline. This is expected to be located in the project's  ||
## directory. The full path to this file should be specified  ||
## in the main config file in the "CONFIG" variable.          ||
##                                                            ||
## The project's directory should be specified in the "dir"   ||
## variable.                                                  ||
##                                                            ||
## ===========================================================##

## Specify the full path to the project's directory ##
dir <- "/path/to/project/directory"

metaDir <- paste0(dir, "/0_metadata")
dataDir <- paste0(dir, "/1_raw")
fastQCDir <- paste0(dir, "/1_raw/fastqc")
trimDir <- paste0(dir, "/2_trimmed")
alignedDir <- paste0(dir, "/3_aligned")
peakDir <- paste0(dir, "/4_calledPeaks")
qcDir <- paste0(alignedDir, "/QCOutput")
sampleSheet <- paste0(metaDir, "/sampleSheet.csv")
analysisDir <- paste0(dir, "/6_analysis")


## create colourblind friendly palette
colorBlindGrey8 <- c(
    "#009E73", "#CC79A7", "#D55E00", "#999999",
    "#F0E442", "#0072B2", "#E69F00", "#56B4E9"
)

## Cell types (list of cell types ordered alphabetically
cellTypes <- c("IRF8", "NEUN", "SOX10", "TN")

## Thresholds applied in the first stage of QC
nrfThres <- 0.7
nucfThres <- 0.3
monoThres <- 0.2
dipPThres <- 0.5
alignThres <- 0.8
uniqReads <- 2
pbc1Thres <- 0.7
pbc2Thres <- 1
peakThres <- 0
fripThres <- 0
