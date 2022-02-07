## this script generates the fragment distribution and calculates summary statistics to assess the periodicity
## given the heavy computational burden, QC is split into batches of 10 samples 

## EXECUTION
# Rscript ./ATACSeq/preprocessing/4_fragmentDistribution.sh <aligned folder> <batch number>
# where 
# <aligned folder> is the folder contained aligned, filtered bam files
# <batch number> is the batch of samples to process
# script needs to be executed from <git repo>/sequencing/


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

args <- commandArgs()
alignedDir<-"/gpfs/mrc0/projects/Research_Project-MRC190311/ATACSeq/rizzardi/3_aligned/" #args[6]
batchNum<- 0 #as.numeric(args[7]) ## nb starts from 0

library(ATACseqQC)
library(diptest)
library(ptest)

## get filepaths of aligned QC'd bam file
aQCFiles<-list.files(alignedDir, pattern = "_depDup_q30.bam$", recursive = TRUE, full.names = TRUE)
aQCSampleNames<-gsub("_depDup_q30.bam", "", basename(aQCFiles))

## filter to subset of samples
index<-c(1:10)+(batchNum*10)
## if number of samples is not a function of ten adjust index
index<-index[index %in% 1:length(aQCFiles)]
nSamples <- length(index)
	
if(length(aQCFiles) > 0){

	## create summary of fragment size using filtered aligned files
	fragSizeHist<-fragSizeDist(aQCFiles, aQCSampleNames)
		
	fragSizeNorm <-lapply(fragSizeHist,standardizeValues)
	## convert to ratios of nucleosomefree, mono, bi, tri etc
	propNucleosomes<-lapply(fragSizeNorm,calcNucleoProps)
	propNucleosomes<-matrix(data = unlist(propNucleosomes), ncol = 5, byrow = TRUE)
	rownames(propNucleosomes)<-aQCSampleNames
	
	## test for multimodality
	diptestStats<-cbind(unlist(lapply(fragSizeNorm, function(y) { dip.test(y)$statistic })), unlist(lapply(fragSizeNorm, function(y) { dip.test(y)$p.value })))
	colnames(diptestStats)<-c("D", "p.value")

	## test for periodicity
	periodTestStats <- cbind(unlist(lapply(fragSizeNorm, function(y) { ptestg(y,method="Fisher")$obsStat })),unlist(lapply(fragSizeNorm, function(y) { ptestg(y,method="Fisher")$pvalue })), unlist(lapply(fragSizeNorm, function(y) { ptestg(y,method="Fisher")$freq })))
	colnames(periodTestStats)<-c("obsStat", "p.value", "freq")
	

	save(fragSizeHist, propNucleosomes, diptestStats, periodTestStats, file = paste0(alignedDir, "/QCOutput/FragmentDistribution_Batch", batchNum, ".rdata"))
}
