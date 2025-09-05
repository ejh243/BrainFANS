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
dir<- "/path/to/project/directory"

metaDir<- paste0(dir, "/0_metadata")
dataDir<-paste0(dir, "/1_raw")
fastQCDir<-paste0(dir, "/1_raw/fastqc")
trimDir<-paste0(dir, "/2_trimmed") 
alignedDir<-paste0(dir, "/3_aligned")
peakDir<- paste0(dir, "/4_calledPeaks")
countsDir <- paste0(dir, "/5_countPeaks")
qcDir<-paste0(alignedDir, "/QCOutput" )
analysisDir<-paste0(dir, "/6_analysis")
sexChrCountsDir <-paste0(dir,"/5_countPeaks/sexChrCountsInPeaks") #(change to user's input)
samplesPeaks <- paste0(peakDir,"/path/to/sample/peaks/")
groupsPeaks <- paste0(peakDir, "/path/to/group/peaks/")
groupsCounts <- paste0(countsDir, "/path/to/group/counts/")
ccPeaks <- paste0(peakDir, "/path/to/case-control/peaks/")
ccCounts <- paste0(peakDir, "/path/to/case-control/counts/")


## Samples information files (change to user's input) ##
sampleSheet<-paste0(metaDir, "/sampleSheet.csv")
samplesList <- paste0(metaDir, "/samples.txt")

## QC stages outputs ##
# Stage 1 results files
qc1stats <- paste0(metaDir, "/qc1Stats.csv")
qc1pass <- paste0(metaDir, "/passStage1Status.csv")
qc1passList <- paste0(metaDir, "/passStage1SampleList.txt")

# Stage 2 results files
qc2pass <- paste0(metaDir, "/passStage2Status.csv")
qc2passList <- paste0(metaDir, "/passStage2SampleList.txt")

qc12stats <- paste0(metaDir, "/passS1S2Status.csv")
groupAnalysisSamples <- paste0(metaDir, "/samplesGroupAnalysis.csv")

# Stage 3 results files
qc3pass <- paste0(metaDir, "/passStage3Status.csv")
qc3passList <- paste0(metaDir, "/passStage3SampleList.txt")
qc3pass_corrected <- paste0(metaDir, "/passStage3Status.csv")

caseControlSamples <- paste0(metaDir, "/samplesCaseControl.csv")

## create colourblind friendly palette
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")
                       
## Thresholds applied in the first stage of QC (change to user's input)
nrfThres <- 0.7
nucfThres <- 0.3
monoThres <- 0.2
dipPThres <-0.5
alignThres <- 0.8
uniqReads <- 2
pbc1Thres <- 0.7
pbc2Thres <- 1
peakThres <- 0
fripThres <- 0

## Peak calling thresholds (change to user's input)
macs_thres_group <- "1e-2"
macs_thres_sample <- "5e-2"

## Cell-type check thresholds (change to user's input)
## Number of top differential peaks to do cell-type check
diffN <- 5000
## Silhouette score threshold to mark failed samples in cell-type check
silScore <- -0.1
                 
# Cell types (list of cell types ordered alphabetically (change to user's input)
cellTypes <- c("IRF8", "NEUN", "SOX10", "TN")

