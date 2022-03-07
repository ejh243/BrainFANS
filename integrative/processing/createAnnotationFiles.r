## need BEDTools and BEDOPs available on the command line

setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/")

library(data.table)
library(bedr)
library(dplyr)


## load ATAC peakSets
allpeaks<-read.table("ATACSeq/CalledPeaks/AllData/OCT2020/AllFractions_peaks.broadPeak", stringsAsFactors = FALSE)
allpeaks<-allpeaks[order(allpeaks[,1], allpeaks[,2]),]
allRegions.atac<-bed2index(allpeaks[,1:3])

totalpeaks<-read.table("ATACSeq/CalledPeaks/AllData/OCT2020/T_peaks.broadPeak", stringsAsFactors = FALSE)
totalpeaks<-totalpeaks[order(totalpeaks[,1], totalpeaks[,2]),]
totalRegions.atac<-bed2index(totalpeaks[,1:3])

neunpeaks<-read.table("ATACSeq/CalledPeaks/AllData/OCT2020/N_peaks.broadPeak", stringsAsFactors = FALSE)
neunpeaks<-neunpeaks[order(neunpeaks[,1], neunpeaks[,2]),]
neunRegions.atac<-bed2index(neunpeaks[,1:3])

soxpeaks<-read.table("ATACSeq/CalledPeaks/AllData/OCT2020/S_peaks.broadPeak", stringsAsFactors = FALSE)
soxpeaks<-soxpeaks[order(soxpeaks[,1], soxpeaks[,2]),]
soxRegions.atac<-bed2index(soxpeaks[,1:3])


### load 5hmC peaks

allpeaks<-read.table("CEGX/calledPeaks/All_peaks.narrowPeak", stringsAsFactors = FALSE)
allpeaks<-allpeaks[order(allpeaks[,1], allpeaks[,2]),]
allRegions.5hmc<-bed2index(allpeaks[,1:3])

totalpeaks<-read.table("CEGX/calledPeaks/Total_peaks.narrowPeak", stringsAsFactors = FALSE)
totalpeaks<-totalpeaks[order(totalpeaks[,1], totalpeaks[,2]),]
totalRegions.5hmc<-bed2index(totalpeaks[,1:3])

neunpeaks<-read.table("CEGX/calledPeaks/NeuN_peaks.narrowPeak", stringsAsFactors = FALSE)
neunpeaks<-neunpeaks[order(neunpeaks[,1], neunpeaks[,2]),]
neunRegions.5hmc<-bed2index(neunpeaks[,1:3])

soxpeaks<-read.table("CEGX/calledPeaks/Sox10_peaks.narrowPeak", stringsAsFactors = FALSE)
soxpeaks<-soxpeaks[order(soxpeaks[,1], soxpeaks[,2]),]
soxRegions.5hmc<-bed2index(soxpeaks[,1:3])

dnegpeaks<-read.table("CEGX/calledPeaks/DoubleNegative_peaks.narrowPeak", stringsAsFactors = FALSE)
dnegpeaks<-dnegpeaks[order(dnegpeaks[,1], dnegpeaks[,2]),]
dnegRegions.5hmc<-bed2index(dnegpeaks[,1:3])


## create a background list of all peaks across all peak sets
allReg<-bedr.merge.region(bedr.sort.region(c(allRegions.atac, totalRegions.atac, neunRegions.atac, soxRegions.atac, allRegions.5hmc, totalRegions.5hmc, neunRegions.5hmc, soxRegions.5hmc,dnegRegions.5hmc)))

for(chr in 1:22){
	annot<-read.table(gzfile(paste("/gpfs/mrc0/projects/Research_Project-MRC190311/References/LDScore/resources/grch38/baselineLD_v2.2/baselineLD.", chr,".annot.gz", sep = "")),header = TRUE)
	snps<-fread(paste("/gpfs/mrc0/projects/Research_Project-MRC190311/References/LDScore/resources/grch38/plink_files/1000G.EUR.hg38.", chr,".bim", sep = ""))

	## to create a single annotation of SNPs in baseline annotation filter to binary annotations
	annot<-select_if(annot, is.integer)
	
	## extract SNP info
	snpBed<-annot[,c("CHR", "BP", "BP")]
	snpBed<-snpBed[order(snpBed$CHR, snpBed$BP),]
	colnames(snpBed)<-c("chr", "start", "end")
	snpBed$end<-snpBed$end+1
	snpBed$chr<-paste("chr", snpBed$chr, sep = "")
	class(snpBed$end)<-"integer"
	snpRegions<-bed2index(as.data.frame(snpBed))
	
	## classify SNPs in atac peaks
	atacALL<-in.region(snpRegions, allRegions.atac)
	atacTOTAL<-in.region(snpRegions, totalRegions.atac)
	atacNEUN<-in.region(snpRegions, neunRegions.atac)
	atacSOX<-in.region(snpRegions, soxRegions.atac)
	hmcALL<-in.region(snpRegions, allRegions.5hmc)
	hmcTOTAL<-in.region(snpRegions, totalRegions.5hmc)
	hmcNEUN<-in.region(snpRegions, neunRegions.5hmc)
	hmcSOX<-in.region(snpRegions, soxRegions.5hmc)
	hmcDNEG<-in.region(snpRegions, dnegRegions.5hmc)
	overlapany<-in.region(snpRegions, allReg)
	mergeOverlaps<-cbind(atacALL,atacTOTAL,atacNEUN, atacSOX, hmcALL, hmcTOTAL, hmcNEUN, hmcSOX, hmcDNEG)
	## createcolumn of overlap with other regulatory features but not in any atac peak
	outsidepeaks<-as.numeric(rowSums(annot[,-c(1:3)]) > 0)
		
	annOut<-data.frame(snps$V1, snps$V4, snps$V2, snps$V3, mergeOverlaps, outsidepeaks, annot$base)
	colnames(annOut)<-c("CHR", "BP", "SNP", "CM","All_ATAC", "Total_ATAC", "NeuN_ATAC", "Sox10_ATAC", "All_5hmc", "Total_5hmc", "NeuN_5hmc", "Sox10_5hmc","DoubleNeg_5hmc", "AnyRegFeature", "AllSNPs")
	mode(annOut$All_ATAC)<-"numeric"
	mode(annOut$Total_ATAC)<-"numeric"
	mode(annOut$NeuN_ATAC)<-"numeric"
	mode(annOut$Sox10_ATAC)<-"numeric"
	mode(annOut$All_5hmc)<-"numeric"
	mode(annOut$Total_5hmc)<-"numeric"
	mode(annOut$NeuN_5hmc)<-"numeric"
	mode(annOut$Sox10_5hmc)<-"numeric"
	mode(annOut$DoubleNeg_5hmc)<-"numeric"

	write.table(annOut, paste("LDScoreRegression/NeuralAnnotations/NeuralCellRegulatoryPeaks.", chr, ".annot", sep = ""), quote = FALSE, row.names = FALSE)
}