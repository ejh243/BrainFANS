## create summary of cohort that passed DNAm QC

args<-commandArgs(trailingOnly = TRUE)
source(args[1])

library(bigmelon)

setwd(dataDir)

load(normData)

## filter samples
QCSum<-read.csv(paste0(qcOutFolder,"PassQCStatusAllSamples.csv"), stringsAsFactors = FALSE, row.names = 1)
passQC<-rownames(QCSum)[which(QCSum$predLabelledCellType == "TRUE")]
pheno<-pheno[match(passQC, pheno$Basename),]

cellTypes<-unique(pheno$Cell.type)
cellTypes<-cellTypes[!is.na(cellTypes)]
## sort so colours lines up correctly
cellTypes<-sort(cellTypes)
cellCols<-c("darkgreen", "darkblue", "darkmagenta", "deeppink", "darkgray") ## assumes celltypes are order alphabetically
names(cellCols)<-cellTypes

## create dataframe with unique entry for each individual
uniqueIDs<-unique(pheno[,c("Indidivual.ID", "Sex", "Age", "Phenotype", "Tissue.Centre")])

## age distributions
par(mfrow = c(2,3))
age_lim<-range(uniqueIDs$Age,na.rm = TRUE)
hist_breaks<-seq(age_lim[1],age_lim[2], length.out = 25)
hist(pheno$Age, hist_breaks, ylab = "nIndividuals",xlab = "Age",main = "All")
for(each in cellTypes){
	index<-which(pheno$Cell.type == each)
	hist(pheno$Age[index], hist_breaks, col = cellCols[each], ylab = "nIndividuals",xlab = "Age",main = each)
}

## sex distribution
table(pheno$Cell.type, pheno$Sex)

## brain bank distribution
table(pheno$Cell.type, pheno$Tissue.Centre)

## schizophrenia status distribution
table(pheno$Cell.type, pheno$Phenotype)

## compare age prediction by cell type
par(mfrow = c(2,3))
for(each in cellTypes){
	index<-which(pheno$Cell.type == each)
	model<-lm(pheno$Age~pheno$DNAmAge, subset = index)
	plot(pheno$DNAmAge[index], pheno$Age[index], xlab = "Horvath Clock", ylab = "Reported", main=each, pch=16, col=cellCols[each])
	title(main = paste("r = ", signif(cor(pheno$DNAmAge[index], pheno$Age[index], use = "pairwise.complete.obs"),3)), line = 0.5, adj = 1)
	title(main = paste("RMSE = ", signif(sqrt(median(abs(pheno$DNAmAge[index]-pheno$Age[index])^2, na.rm = TRUE)),3)), line = 2, adj = 1)
	abline(model, lty = 2)
	abline(a = 0, b = 1)
}

for(each in cellTypes){
	index<-which(pheno$Cell.type == each)

	model<-lm(pheno$Age~pheno$CCDNAmAge, subset = index)
	plot(pheno$CCDNAmAge[index], pheno$Age[index], xlab = "Cortical Clock", ylab = "Reported", main=each, pch=16, col=cellCols[each])
	title(main = paste("r = ", signif(cor(pheno$CCDNAmAge[index], pheno$Age[index], use = "pairwise.complete.obs"),3)), line = 0.5, adj = 1)
	title(main = paste("RMSE = ", signif(sqrt(median(abs(pheno$CCDNAmAge[index]-pheno$Age[index])^2, na.rm = TRUE)),3)), line = 2, adj = 1)
	abline(model, lty = 2)
	abline(a = 0, b = 1)
}

