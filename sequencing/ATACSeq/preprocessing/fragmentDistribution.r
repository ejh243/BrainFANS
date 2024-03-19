## =========================================================================================================================##
##                             ATAC-seq pipeline STEP 2.1: fragment distribution and post-alignment metrics                 ##
## =========================================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/preprocessing/fragmentDistribution.sh <project> <batch number>                   ||
## - execute from scripts directory                                                                                         ||
##                                                                                                                          ||
## DESCRIPTION: This script uses the ATACseqQC R package to generate the fragment distribution and calculate some           ||
## summary statistics to assess the periodicity                                                                             ||
##                                                                                                                          ||
## INPUTS:                                                                                                                  ||
## - <project> : project on which analysis is being run                                                                    ||
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
  list_plots <- vector("list", 10)
  samples <- ls(listSamples)
  for(i in 1:length(listSamples)) {
    df <-data.frame(listSamples[i])
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

args <- commandArgs()
project <-args[6]
source(paste0("/lustre/projects/Research_Project-MRC190311/ATACSeq/",project,"/config.r")
batchNum<-as.numeric(args[7]) ## nb starts from 0

library(ATACseqQC)
library(diptest)
library(ptest)
library(plyr)
library(ggplot2)
library(ggpubr)

## get filepaths of aligned indexed QC'd bam file
aQCFiles<-list.files(alignedDir, pattern = ".filt.nodup.bam.bai$", recursive = TRUE, full.names = TRUE)
aQCFiles<-gsub(".bai", "", aQCFiles)
aQCSampleNames<-gsub(".filt.nodup.bam", "", basename(aQCFiles))

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
	save(fragSizeHist, propNucleosomes, diptestStats, periodTestStats, file = paste0(alignedDir, "/QCOutput/FragmentDistribution_Batch", batchNum, ".rdata"))
 
  ## output Fragment size distribution plots for the corresponding batch of samples
  
  ## split samples in two groups for plotting purposes
  
  if(length(fragSizeHist) < 10){
    fragSizeHist.1<-fragSizeHist[c(1:length(fragSizeHist))]
    fragSizeHist.2<-fragSizeHist.1
    
  } else {
    fragSizeHist.1<-fragSizeHist[c(1:5)]
    fragSizeHist.2<-fragSizeHist[c(6:10)]
    
  }
  
  BatchSamples.1 <-fdsPlots(fragSizeHist.1)
  BatchSamples.2 <-fdsPlots(fragSizeHist.2)

  axis1 <- seq(from=0, to=max(BatchSamples.1[,1]), by=100)
  p.01 <- ggplot(BatchSamples.1, aes(x=BatchSamples.1[,1], y=BatchSamples.1[,2], group=sample, col=sample, fill=sample)) + geom_line()+
  labs(x="Fragment size (bp)",y="Frequency", title=paste0("Fragment size distribution batch ",batchNum," ,set 1")) +scale_x_continuous(breaks=axis1)
  
  axis2 <- seq(from=0, to=max(BatchSamples.2[,1]), by=100)
  p.02 <- ggplot(BatchSamples.2, aes(x=BatchSamples.2[,1], y=BatchSamples.2[,2], group=sample, col=sample, fill=sample)) + geom_line()+scale_x_continuous(breaks=axis2)+
  labs(x="Fragment size (bp)",y="Frequency", title=paste0("Fragment size distribution batch ",batchNum, " ,set 2"))
  
  ## Plots for batch will be output in a single pdf
  plots<-ggarrange(p.01, p.02,  
          ncol = 1, nrow = 2)

  pdf(file.path(paste0(alignedDir, "/QCOutput/FSD_batch", batchNum, ".pdf")), width = 9, height = 7)
  print(plots)
  dev.off()
  
}


