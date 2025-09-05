## =========================================================================================================================##
##                             ATAC-seq pipeline STEP 2.1: fragment distribution and post-alignment metrics                 ##
## =========================================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/fragmentDistribution.r <project> <batch number>                         ||
## - execute from scripts directory                                                                                         ||
##                                                                                                                          ||
## DESCRIPTION: This script uses the ATACseqQC R package to generate the fragment distribution and calculate some           ||
## summary statistics to assess the periodicity                                                                             ||
##                                                                                                                          ||
## INPUTS:                                                                                                                  ||
## - <project> : project on which analysis is being run                                                                     ||
## - <batch number>: number of samples to process at a time, should be <10                                                  ||
##                                                                                                                          ||
## OUTPUTS:                                                                                                                 ||
## - FragmentDistribution_Batch.rdata, FSD_batch<batchNumber>.pdf                                                           ||
##                                                                                                                          ||
## REQUIRES:                                                                                                                ||
## - R/4.2.1-foss-2022a, libraries: ATACseqQC, diptest, ptest, plyr ,ggplot2                                                ||
## =========================================================================================================================##

## ==========##
## FUNCTIONS ##
## ==========##

standardizeValues<-function(frag.len){
    x <- 1:1010
    frag.len <- frag.len[match(x, names(frag.len))]
    frag.len[is.na(frag.len)] <- 0
    y <- frag.len / sum(frag.len)
    y <- as.numeric(y)
    names(y) <- x
	return(y)
	}
	
calcNucleoProps<-function(fraglen.stand){
	freeSum<-sum(fraglen.stand[1:150])
	monoSum<-sum(fraglen.stand[151:300])
	diSum<-sum(fraglen.stand[301:450])
	triSum<-sum(fraglen.stand[451:600])
	otherSum<-sum(fraglen.stand[-c(1:600)])
	return(c(freeSum, monoSum, diSum, triSum, otherSum))
}

fdsPlots <- function(listSamples){
  counter = 1
  list_plots <-list()
  samples <- names(listSamples)
  for(i in 1:length(listSamples)) {
    df <-data.frame(listSamples[[i]])
    colnames(df)<- c("Var1", "Freq")
    df$sample <-samples[i]
    list_plots[[i]] <- df
    counter = counter + 1
  }
  bigdf <-ldply(list_plots, data.frame)
  bigdf[,1]<- as.integer(as.character(bigdf[,1]))
  return(bigdf)
}

## ==========##
##   SET-UP  ##
## ==========##

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)
batchNum<-as.numeric(args[2]) ## nb starts from 0

options(scipen=5)
suppressWarnings(suppressPackageStartupMessages({
  library(ATACseqQC)
  library(diptest)
  library(ptest)
  library(plyr)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(dplyr)
}))

## get filepaths of aligned indexed QC'd bam file
samples<-read.table(file.path(samplesList))[,1]
aQCFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam.bai$", recursive = TRUE, full.names = TRUE)
aQCFiles<-gsub(".bai", "", aQCFiles)
aQCFiles <- aQCFiles[match(samples,gsub(paste0(alignedDir,"/"),"",gsub(".filt.nodup.bam", "", aQCFiles)))]
aQCSampleNames<-gsub(".filt.nodup.bam", "", basename(aQCFiles))
aQCSampleNames <- na.omit(aQCSampleNames)
aQCFiles <- na.omit(aQCFiles)
## filter to subset of samples specified by array number
index<-c(1:10)+(batchNum*10)
## if number of samples is not a function of ten adjust index
index<-index[index %in% 1:length(aQCFiles)]
nSamples <- length(index)


## =============##
##   CALCULATE  ##
## =============##

if(nSamples > 0){

	## create summary of fragment size using filtered aligned files
  print("Samples in batch: ")
  print(aQCSampleNames[index])
	fragSizeHist<-fragSizeDist(aQCFiles[index], aQCSampleNames[index])
		
	fragSizeNorm <-lapply(fragSizeHist,standardizeValues)
 
	## convert to ratios of nucleosomefree, mono, bi, tri etc
	propNucleosomes<-lapply(fragSizeNorm,calcNucleoProps)
	propNucleosomes<-matrix(data = unlist(propNucleosomes), ncol = 5, byrow = TRUE)
	rownames(propNucleosomes)<-aQCSampleNames[index]
	
	## test for multimodality
	diptestStats<-cbind(unlist(lapply(fragSizeNorm, function(y) { dip.test(y)$statistic })), unlist(lapply(fragSizeNorm, function(y) { dip.test(y)$p.value })))
	colnames(diptestStats)<-c("D", "p.value")

	## test for periodicity
	periodTestStats <- cbind(unlist(lapply(fragSizeNorm, function(y) { ptestg(y,method="Fisher")$obsStat })),unlist(lapply(fragSizeNorm, function(y) { ptestg(y,method="Fisher")$pvalue })), unlist(lapply(fragSizeNorm, function(y) { ptestg(y,method="Fisher")$freq })))
	colnames(periodTestStats)<-c("obsStat", "p.value", "freq")
	
  ## Results are saved in a rdata file 
	save(fragSizeHist, propNucleosomes, diptestStats, periodTestStats, file = paste0(qcDir, "/FragmentDistribution_Batch_", batchNum, ".rdata"))
 
  ## output Fragment size distribution plots for the corresponding batch of samples
  
  ## split samples in two groups for plotting purposes
  BatchSamples <- fdsPlots(fragSizeHist)
  BatchSamples.split <- split(BatchSamples,BatchSamples$sample)
  p<- NULL
  
  for(i in 1:length(BatchSamples.split)){
    fdValues <- BatchSamples.split[[i]]
    axis1 <- seq(from=0, to=max(fdValues[,1]), by=100)
    p[[i]] <- ggplot(fdValues, aes(x=Var1, y=Freq)) + geom_line()+theme_bw()+
  labs(x="Fragment size (bp)",y="Frequency", title=names(BatchSamples.split)[i]) +scale_x_continuous(breaks=axis1, limits = c(0, 800))
  }

  pdf(file.path(paste0(qcDir, "/FSD_batch_", batchNum, ".pdf")), width = 10, height = 5)
  print(p)
  dev.off()
  
  
}


