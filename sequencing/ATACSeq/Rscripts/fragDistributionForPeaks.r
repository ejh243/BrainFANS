## =========================================================================================================================##
##                   ATAC-seq pipeline STEP 7.1: Fragment size distribution of samples for group peak calling               ##
## =========================================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/fragDistributionForPeaks.r <project> <cell-group>                       ||
## - execute from scripts directory                                                                                         ||
##                                                                                                                          ||
## DESCRIPTION: This script uses the ATACseqQC R package to generate the fragment distribution and calculate some           ||
## summary statistics to assess the periodicity                                                                             ||
##                                                                                                                          ||
## INPUTS:                                                                                                                  ||
## - <project> : project on which analysis is being run                                                                     ||
## - <cell-group>: number of samples to process at a time, should be <10                                                    ||
##                                                                                                                          ||
## OUTPUTS:                                                                                                                 ||
## - FragmentDistribution_<cell-group>.pdf in 5_countPeaks folder                                                           ||
##                                                                                                                          ||
## REQUIRES:                                                                                                                ||
## - R/4.2.1-foss-2022a, libraries: ATACseqQC,plyr ,ggplot2                                                                 ||
## =========================================================================================================================##

## ==========##
## FUNCTIONS ##
## ==========##
options(scipen=5)

library(ggplot2)
library(ATACseqQC)
library(plyr)

args <- commandArgs()
configR <-source(args[6])
cf <- args[7]

## ==========##
##   SET-UP  ##
## ==========##

fdsPlots <- function(listSamples){
  list_plots <-list()
  samples <- names(listSamples)
  for(i in 1:length(listSamples)) {
    df <-data.frame(listSamples[i])
    colnames(df)<- c("Var1", "Freq")
    df$sample <-samples[i]
    list_plots[[i]] <- df
  }
  bigdf <-ldply(list_plots, data.frame)
  bigdf[,1]<- as.integer(as.character(bigdf[,1]))
  return(bigdf)
}

samples <- read.table(file.path(paste0(metaDir, "/samplesForGroupAnalysisOrdered_",cf,".txt")), stringsAsFactors = FALSE)$V1
aQCFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam.bai$", recursive = TRUE, full.names = TRUE)
aQCFiles<-gsub(".bai", "", aQCFiles)
aQCFiles <- aQCFiles[match(samples,gsub(paste0(alignedDir,"/"),"",gsub(".filt.nodup.bam", "", aQCFiles)))]
aQCSampleNames<-gsub(".filt.nodup.bam", "", basename(aQCFiles))
aQCSampleNames <- na.omit(aQCSampleNames)
aQCFiles <- na.omit(aQCFiles)
nSamples <- length(samples)

if(nSamples > 0){

  print("Samples for fragment size distribution: ")
  print(samples)
	fragSizeHist<-fragSizeDist(aQCFiles, samples)
 
  BatchSamples <- fdsPlots(fragSizeHist)
  BatchSamples.split <- split(BatchSamples,BatchSamples$sample)
  p<- NULL
  
  for(i in 1:length(BatchSamples.split)){
    fdValues <- BatchSamples.split[[i]]
    axis1 <- seq(from=0, to=max(fdValues[,1]), by=100)
    p[[i]] <- ggplot(fdValues, aes(x=Var1, y=Freq)) + geom_line()+theme_bw()+
  labs(x="Fragment size (bp)",y="Frequency", title=names(BatchSamples.split)[i]) +scale_x_continuous(breaks=axis1)
  }
  
  pdf(file.path(paste0(dir,"/5_countPeaks/FragmentDistribution_", cf, ".pdf")), width = 10, height = 5)
  print(p)
  dev.off()
  
}