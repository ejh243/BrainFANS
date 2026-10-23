## ===========================================================##
##                   PIPPA: Example config R file             ##
## ===========================================================##
##                                                            ||
## This is an example config R file for PIPPA. This is        ||
## expected to be located in the project's  directory. The    ||
## full path to this file should be specified in the main     ||
## in the main config file in the "CONFIGR" variable.         ||
##                                                            ||
## The project's directory should be specified as "mainDir"   ||
## variable.                                                  ||
##                                                            ||
## ===========================================================##

## Specify path to pipeline's config directory ##
configDir <- "./config/"

## Specify the full path to the project's directory ##
mainDir<- "/path/to/project/directory"

metaDir<- paste0(mainDir, "/0_metadata")
dataDir<-paste0(mainDir, "/1_raw")
fastQCDir<-paste0(mainDir, "/1_raw/fastqc")
trimDir<-paste0(mainDir, "/2_trimmed") 
alignedDir<-paste0(mainDir, "/3_aligned")
peakDir<- paste0(mainDir, "/4_calledPeaks")
countsDir <- paste0(mainDir, "/5_countPeaks")
qcDir<-paste0(alignedDir, "/QCOutput")
analysisDir<-paste0(mainDir, "/6_analysis")
sexChrCountsDir <-paste0(mainDir,"/5_countPeaks/sexChr/") 
genoFolder<-paste0(alignedDir, "/genotypeConcordance/")

## Peaks directories (change to user's input) ##
samplesPeaks <- paste0(peakDir,"/samplePeaks/")
groupsPeaks <- paste0(peakDir, "/groupCheck/")
peaksAnalysis <- paste0(peakDir, "/peaksCuration/peaksAnalysis/")
peaksCutoff <- paste0(peakDir,"/peaksCuration/cutoff/")
peaksCurated <-  paste0(peakDir,"/peaksCuration/curatedPeaks/")

## Counts directories (change to user's input) ##
groupsCounts <- paste0(countsDir, "/checkCTCounts/")
downstreamCounts <- paste0(countsDir, "/caseControl/")

## create colourblind friendly palette
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")
##change as desired
col_vals_plots <- c("SOX10+" = "#6B8E23","TN" = "#C77CFF","NeuN+" = "hotpink","IRF8+" = "orange")
                       
## Thresholds applied in the first stage of QC (change to user's input)
nrfThres <- 0.7
nucfThres <- 0.25
monoThres <- 0.2
alignThres <- 0.8
uniqReads <- 50
pbc1Thres <- 0.7
pbc2Thres <- 1

## Peak calling thresholds (change to user's input)
macs_thres_group <- "5e-2"
macs_thres_sample <- "5e-2"

## Cell-type check thresholds (change to user's input)
## Number of top differential peaks to do cell-type check
diffN <- 1000

## Silhouette score threshold to mark failed samples in cell-type check
silScore <- 0
                 
# Cell types (list of cell types ordered alphabetically (change to user's input)
cellTypes <- c("IRF8+", "NeuN+", "SOX10+", "TN")

## clinical groups (optional)
clinicalGroups <- c("case", "control")
## brain banks groups (optional)
brainBanks <- c("A", "B","C")

## Samples information files (change to user's input) ##
sampleSheet<-paste0(metaDir, "/sampleSheet.csv")
samplesList <- paste0(metaDir, "/samples.txt")

## Files for QC (change to user's input)##
fripFile <- paste0(peakDir, "/QCOutput/FRIP_all_samples.csv")
tssEnrich <- paste0(qcDir, "/TSS_enrichment.tsv")
processSumFileS1 <- paste0(metaDir, "/summariseSampleProcessingProgress.csv")
fastqcFile <- paste0(fastQCDir, "/multiqc/multiqc_data/multiqc_fastqc.txt")
readCountsFile <- paste0(qcDir,"/countMTReads.txt")
alignmentFile <- paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt")

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

# Overall QC results
allQCresults <- paste0(metaDir, "/passAllQCStatus.csv")
allQCstats <- paste0(metaDir, "/resultsAllQC.csv")
downAnalysisSamples <- paste0(metaDir, "/samplesDownAnalysis.csv")

## Counts in peaks files for STEP 8 (change to user's input) ##
groupCountsAll <- paste0(groupsCounts,"/peakCounts_all_",macs_thres_group,"_passS12.narrowPeak.rdata")
diffPeaksCountsAll <- paste0(groupsCounts, "/cellType_diffPeaks_edgeR_",macs_thres_group,"_passS12.narrowPeak.rdata")

## Counts in curated peaks for STEP 10 (change to user's input) ##
groupCountsClean <- paste0(groupsCounts,"/peakCounts_all_",macs_thres_group,"_passAll.narrowPeak.rdata")
diffPeaksCountsClean <- paste0(groupsCounts, "/cellType_diffPeaks_edgeR_",macs_thres_group,"_passAll.narrowPeak.rdata")
