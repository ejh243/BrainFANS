##---------------------------------------------------------------------#
##
## Title: 
##
## Purpose of script: to collate the peak information from the whole-sample peak calling
##
## Author: 
##
## Date Created: 2022-05-27
##
##---------------------------------------------------------------------#
##
## set working directory
setwd("")
## clear the R environment
rm(list=ls()) 

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
project<-"hypothalamus"
source("/lustre/projects/Research_Project-MRC190311/scripts/sequencing/ATACSeq/config/config.r")

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
#good and bad samples
## ALIGNMENT STATISTICS
alignQC<-read.table(file.path(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## FRIP
fripFiles <- list.files(file.path(peakDir, "/QCOutput/subset/"), pattern = "FRIP")


## 6. FRIP
if ( length(fripFiles) != 0 ) {
  fripStats<-read.csv(paste0(peakDir, "/QCOutput/subset/", fripFiles[1]))
  for(each in fripFiles[-1]){
    fripStats<-rbind(fripStats, read.csv(paste0(peakDir, "/QCOutput/subset/", each)))
  }
  
  #fripStats<-fripStats[match(processSum$sampleID, fripStats$SampleName),]
} else { warning('FRIP stats do not appear to have been calculated') }

fripStats$FripSample <- fripStats$readsInSampleTagAlignPeaks/fripStats$tagAlignTotalReads
fripStats$FripSubset <- fripStats$readsInSubsetTagAlignPeaks/fripStats$tagAlignTotalReads

par(mfrow = c(1, 1))
vioplot(fripStats$FripSample, fripStats$FripSubset, ylab = "Fraction reads in peaks", names = c("Sample-specific", "Sample reads in whole peak set"))

vioplot(fripStats$FripSubset[alignQC$overall_alignment_rate > 69], 
        fripStats$FripSubset[alignQC$overall_alignment_rate < 70], ylab = "Fraction reads in peaks",
        names = c("Pass", "Fail"))

plot(fripStats$FripSubset[alignQC$overall_alignment_rate < 70])


sampleSheet
