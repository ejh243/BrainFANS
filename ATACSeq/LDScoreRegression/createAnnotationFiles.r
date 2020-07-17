## need BEDTools and BEDOPs available on the command line
## nb peaks are from hg38 LDscores are hg19

setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/ATACSeq/CalledPeaks/AllData/FANSPaper/")

library(data.table)
library(bedr)

peakFiles<-list.files(".", pattern = "broadPeak")

## load peakSets
allpeaks<-read.table("AllFractions_peaks.broadPeak", stringsAsFactors = FALSE)
allpeaks<-allpeaks[order(allpeaks[,1], allpeaks[,2]),]
allRegions<-bed2index(allpeaks[,1:3])

totalpeaks<-read.table("TOTAL_peaks.broadPeak", stringsAsFactors = FALSE)
totalpeaks<-totalpeaks[order(totalpeaks[,1], totalpeaks[,2]),]
totalRegions<-bed2index(totalpeaks[,1:3])

neunpeaks<-read.table("NEUN_peaks.broadPeak", stringsAsFactors = FALSE)
neunpeaks<-neunpeaks[order(neunpeaks[,1], neunpeaks[,2]),]
neunRegions<-bed2index(neunpeaks[,1:3])

soxpeaks<-read.table("SOX_peaks.broadPeak", stringsAsFactors = FALSE)
soxpeaks<-soxpeaks[order(soxpeaks[,1], soxpeaks[,2]),]
soxRegions<-bed2index(soxpeaks[,1:3])

## create a background list of all peaks across all peak sets
allReg<-bedr.merge.region(bedr.sort.region(c(allRegions, totalRegions, neunRegions, soxRegions)))

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
	overlapALL<-in.region(snpRegions, allRegions)
	overlapTOTAL<-in.region(snpRegions, totalRegions)
	overlapNEUN<-in.region(snpRegions, neunRegions)
	overlapSOX<-in.region(snpRegions, soxRegions)
	overlapany<-in.region(snpRegions, allReg)
	mergeOverlaps<-cbind(overlapALL,overlapTOTAL,overlapNEUN, overlapSOX)
	## createcolumn of overlap with other regulatory features but not in any atac peak
	outsidepeaks<-as.numeric(rowSums(annot[,-c(1:3)]) > 0)
		
	annOut<-data.frame(snps$V1, snps$V4, snps$V2, snps$V3, mergeOverlaps, outsidepeaks, annot$base)
	colnames(annOut)<-c("CHR", "BP", "SNP", "CM","Allpeakset", "Totalpeakset", "NeuNpeakset", "Sox10peakSet", "AnyRegFeature", "AllSNPs")
	mode(annOut$Allpeakset)<-"numeric"
	mode(annOut$Totalpeakset)<-"numeric"
	mode(annOut$NeuNpeakset)<-"numeric"
	mode(annOut$Sox10peakSet)<-"numeric"

	write.table(annOut, paste("LDScoreAnnotation/ATACPeaks.", chr, ".annot", sep = ""), quote = FALSE, row.names = FALSE)
}