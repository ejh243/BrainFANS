## Written by Eilis
## Summarize CNV calls
## identify any overlapping with known SCZ loci or pathogenic CNVs

findOverlapsMinProp<-function(query, subject, pThres){
	# function to find overlap between two Granges based on minimum proportion
	hits <- findOverlaps(query, subject)
	overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
	percentOverlap <- width(overlaps) / width(subject[subjectHits(hits)])
	hits <- hits[percentOverlap > pThres]
	hitsOut <- query[queryHits(hits)]
	mcols(hitsOut)$Locus<-subject$Locus[subjectHits(hits)]
	mcols(hitsOut)$hg38<-subject$hg38[subjectHits(hits)]
	return(hitsOut)
}


library(GenomicRanges)
args<-commandArgs(trailingOnly = TRUE)
fileName<-args[1]
superPop<-args[2]
folder<-dirname(fileName)

setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/SNPdata/CNV/")
dat<-read.table("PennCNVOutput/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv", stringsAsFactors=FALSE)

pheno<-read.table("../MRC2_UpdatePheno.txt")

dat$V2<-as.numeric(gsub("numsnp=", "",dat$V2))
dat$V3<-as.numeric(gsub(",", "", gsub("length=", "",dat$V3)))

table(dat$V4)

## summarise CNV calls

par(mfrow = c(3,2))
hist(dat$V2, xlab = "number of markers", breaks = 35, main = "All")
hist(dat$V3/1000, xlab = "length (kb)", breaks = 35, main = "All")
index.del<-which(dat$V4 == "state2,cn=1")
hist(dat$V2[index.del], xlab = "number of markers", breaks = 35, main = "Deletions")
hist(dat$V3[index.del]/1000, xlab = "length (kb)", breaks = 35, main = "Deletions")
index.dup<-which(dat$V4 == "state5,cn=3")
hist(dat$V2[index.dup], xlab = "number of markers", breaks = 35, main = "Duplications")
hist(dat$V3[index.dup]/1000, xlab = "length (kb)", breaks = 35, main = "Duplications")


## summarise CNVs by person
par(mfrow = c(3,2))
totByPerson<-aggregate(dat$V3/1000, by = list(dat$V5), sum)
muByPerson<-aggregate(dat$V3/1000, by = list(dat$V5), mean)
hist(totByPerson[,2], xlab = "Combined length of CNVs", ylab = "Number of samples", breaks = 35, main = "All")
hist(muByPerson[,2], xlab = "Mean length of CNVs", ylab = "Number of samples", breaks = 35, main = "All")
totByPerson<-aggregate(dat$V3[index.del]/1000, by = list(dat$V5[index.del]), sum)
muByPerson<-aggregate(dat$V3[index.del]/1000, by = list(dat$V5[index.del]), mean)
hist(totByPerson[,2], xlab = "Combined length of CNVs", ylab = "Number of samples", breaks = 35, main = "Deletions")
hist(muByPerson[,2], xlab = "Mean length of CNVs", ylab = "Number of samples", breaks = 35, main = "Deletions")
totByPerson<-aggregate(dat$V3[index.dup]/1000, by = list(dat$V5[index.dup]), sum)
muByPerson<-aggregate(dat$V3[index.dup]/1000, by = list(dat$V5[index.dup]), mean)
hist(totByPerson[,2], xlab = "Combined length of CNVs", ylab = "Number of samples", breaks = 35, main = "Duplications")
hist(muByPerson[,2], xlab = "Mean length of CNVs", ylab = "Number of samples", breaks = 35, main = "Duplications")


## convert to GRanges; do sep for deletions and duplications
index<-which(dat$V4 == "state2,cn=1")
allCNVs.del<-GRanges(dat$V1[index])
mcols(allCNVs.del)$SampleID <- dat$V5[index]
mcols(allCNVs.del)$Type <- dat$V4[index]
index<-which(dat$V4 == "state5,cn=3")
allCNVs.dup<-GRanges(dat$V1[index])
mcols(allCNVs.dup)$SampleID <- dat$V5[index]
mcols(allCNVs.dup)$Type <- dat$V4[index]

## do any overlap known scz cnv loci ## list taken from Rees et al. Br J Psychiatry (merge tables 1 & 2)
sczLoci<-read.csv("../../References/CNV/SCZ_CNVloci.csv", skip = 1, stringsAsFactors = FALSE)
## filter to those significant in MetaAnalysis
sczLoci<-sczLoci[which(sczLoci$significantMeta == "*"),]
sczLoci.hg38.del<-GRanges(sczLoci$hg38[grep("del", sczLoci$Locus)])
mcols(sczLoci.hg38.del)$Locus<-sczLoci$Locus[grep("del", sczLoci$Locus)]
mcols(sczLoci.hg38.del)$hg38<-sczLoci$hg38[grep("del", sczLoci$Locus)]
sczLoci.hg38.dup<-GRanges(sczLoci$hg38[grep("dup", sczLoci$Locus)])
mcols(sczLoci.hg38.dup)$Locus<-sczLoci$Locus[grep("dup", sczLoci$Locus)]
mcols(sczLoci.hg38.dup)$hg38<-sczLoci$hg38[grep("dup", sczLoci$Locus)]

pThres<-0.9 ## set as minimum overlap required
overlapDel<-findOverlapsMinProp(allCNVs.del, sczLoci.hg38.del, pThres)
overlapDup<-findOverlapsMinProp(allCNVs.dup, sczLoci.hg38.dup, pThres)
output<-rbind(data.frame(overlapDel), data.frame(overlapDup))
write.csv(output, "CNVsoverlappingKnownSCZRiskLoci.csv")

## do any overlap known ID cnv loci ## list taken from Rees et al. JAMA Psychiatry 2016 (eTable 2)
idLoci<-read.csv("../../References/CNV/IDCNVLoci.csv", stringsAsFactors = FALSE)
idLoci<-idLoci[which(idLoci$hg38 != ""),] ## 1 region I couldn't lift over
idLoci.hg38.del<-GRanges(idLoci$hg38[grep("del", idLoci$Syndrome)])
mcols(idLoci.hg38.del)$Locus<-idLoci$Syndrome[grep("del", idLoci$Syndrome)]
mcols(idLoci.hg38.del)$hg38<-idLoci$hg38[grep("del", idLoci$Syndrome)]
idLoci.hg38.dup<-GRanges(idLoci$hg38[grep("dup", idLoci$Syndrome)])
mcols(idLoci.hg38.dup)$Locus<-idLoci$Syndrome[grep("dup", idLoci$Syndrome)]
mcols(idLoci.hg38.dup)$hg38<-idLoci$hg38[grep("dup", idLoci$Syndrome)]

pThres<-0.9 ## set as minimum overlap required
overlapDel<-findOverlapsMinProp(allCNVs.del, idLoci.hg38.del, pThres)
overlapDup<-findOverlapsMinProp(allCNVs.dup, idLoci.hg38.dup, pThres)
output<-rbind(data.frame(overlapDel), data.frame(overlapDup))
write.csv(output, "CNVsoverlappingIDRiskLoci.csv")

## do any overlap known pathogenic cnv loci ## list taken from Kendall et al 2017 Biol Psychiatry
pathLoci<-read.table("../../References/CNV/PathogenicCNVLoci.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
pathLoci<-pathLoci[which(pathLoci$hg38 != ""),] ## 1 region I couldn't lift over
pathLoci.hg38.del<-GRanges(pathLoci$hg38[grep("del", pathLoci$CNV.locus)])
mcols(pathLoci.hg38.del)$Locus<-pathLoci$CNV.locus[grep("del", pathLoci$CNV.locus)]
mcols(pathLoci.hg38.del)$hg38<-pathLoci$hg38[grep("del", pathLoci$CNV.locus)]
pathLoci.hg38.dup<-GRanges(pathLoci$hg38[grep("dup", pathLoci$CNV.locus)])
mcols(pathLoci.hg38.dup)$Locus<-pathLoci$CNV.locus[grep("dup", pathLoci$CNV.locus)]
mcols(pathLoci.hg38.dup)$hg38<-pathLoci$hg38[grep("dup", pathLoci$CNV.locus)]

pThres<-0.9 ## set as minimum overlap required
overlapDel<-findOverlapsMinProp(allCNVs.del, pathLoci.hg38.del, pThres)
overlapDup<-findOverlapsMinProp(allCNVs.dup, pathLoci.hg38.dup, pThres)
output<-rbind(data.frame(overlapDel), data.frame(overlapDup))
write.csv(output, "CNVsoverlappingPathogenicCNV.csv")



