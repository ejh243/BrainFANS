## run this script to confirm sorted cell fractions are labelled correctly and data clusters as expected
## use QC metrics for phenotype data!!
## within adult samples do PCA to cluster sample types


library(bigmelon)
library(pheatmap)
library(glmnet)

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
## exclude low intensity; incomplete bisulfite conversion, incorrect sex prediction, discordant with SNP data
load(qcData)
if(sexCheck){
	QCSum<-cbind(QCmetrics$bisulfCon > thresBS,as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex),(QCmetrics$M.median > intenThres & QCmetrics$U.median > intenThres))
	colnames(QCSum)<-c("BSConversion", "SexPrediction", paste0("intens<", intenThres))
} else {
	QCSum<-cbind(QCmetrics$bisulfCon > thresBS,(QCmetrics$M.median > intenThres & QCmetrics$U.median > intenThres))
	colnames(QCSum)<-c("BSConversion",  paste0("intens<", intenThres))
}
if(snpCheck){
	QCSum<-cbind(QCSum, QCmetrics$genoCheck > 0.8)
	colnames(QCSum)[ncol(QCSum)]<-"snpCheck"
}
rownames(QCSum)<-rownames(QCmetrics)
QCSum<-cbind(QCSum, rowSums(QCSum, na.rm = TRUE) == rowSums(!is.na(QCSum)))
colnames(QCSum)[ncol(QCSum)]<-"passS2"

passQC<-rownames(QCSum)[QCSum[,"passS2"]]

sampleSheet<-QCmetrics
sampleSheet<-sampleSheet[match(passQC, sampleSheet$Basename),]

rawbetas<-gfile[,, node = "betas"]
rawbetas<-rawbetas[,match(passQC, colnames(rawbetas))]
auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
rawbetas<-rawbetas[auto.probes,]

cellTypes<-unique(sampleSheet$Cell.type)
cellCols<-rainbow(length(cellTypes))[as.factor(sampleSheet$Cell.type)]

## filter out NAs
rawbetas<-na.omit(rawbetas)

## filter out SNPs
rawbetas<-rawbetas[-grep("rs", rownames(rawbetas)),]

## use Mahalanobis to calculate difference with cell type medians from PCAs

pca <- prcomp(t(rawbetas))
betas.scores = pca$x
colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
betas.pca<-pca$sdev^2/sum(pca$sdev^2)
betas.scores<-betas.scores[,which(betas.pca > 0.01)]

cellMedPCA<-aggregate(betas.scores, by = list(sampleSheet$Cell.type), median, na.rm = TRUE)
rownames(cellMedPCA)<-cellMedPCA[,1]
cellMedPCA<-as.matrix(t(cellMedPCA[,-1]))

## calc cell type covariance matrices
cov_sigmaPCA<-list()
for(each in colnames(cellMedPCA)){
	cov_sigmaPCA[[each]]<-cov(betas.scores[which(sampleSheet$Cell.type == each),], use = "p")
}

mahDistPCA<-matrix(data = 10^9, ncol = ncol(cellMedPCA), nrow = nrow(sampleSheet))
colnames(mahDistPCA)<-colnames(cellMedPCA)
for(each in colnames(cellMedPCA)){
	if(det(cov_sigmaPCA[[each]]) > 1){
		mahDistPCA[,each]<-mahalanobis(betas.scores, cellMedPCA[,each], cov_sigmaPCA[[each]], na.rm = TRUE,tol=1e-20)
	}
}

closestCellTypePCA<-colnames(mahDistPCA)[unlist(apply(mahDistPCA, 1, which.min))]
predLabelledCellType<-sampleSheet$Cell.type == closestCellTypePCA

sampleSheet<-cbind(sampleSheet, closestCellTypePCA, predLabelledCellType)
write.csv(sampleSheet[which(predLabelledCellType == "FALSE"),], paste0(qcOutFolder, "SamplesPredictedDiffCellTypePCAMahDist.csv"))


## instead calculate outliers

cellMeanPCA<-aggregate(betas.scores, by = list(sampleSheet$Cell.type), mean, na.rm = TRUE)
rownames(cellMeanPCA)<-cellMeanPCA[,1]
cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))

cellSDPCA<-aggregate(betas.scores, by = list(sampleSheet$Cell.type), sd, na.rm = TRUE)
rownames(cellSDPCA)<-cellSDPCA[,1]
cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))

lowerBound<-cellMeanPCA-2*cellSDPCA
upperBound<-cellMeanPCA+2*cellSDPCA



pcaCellClassify<-rep(NA, nrow(sampleSheet))
for(i in 1:length(cellTypes)){
	pcaCellClassify[which(betas.scores[,1] < upperBound[1,i] & betas.scores[,1] > lowerBound[1,i] & betas.scores[,2] < upperBound[2,i] & betas.scores[,2] > lowerBound[2,i])]<-colnames(lowerBound)[i]
}

sampleSheet<-cbind(sampleSheet, pcaCellClassify)

write.csv(sampleSheet[which(pcaCellClassify !=  sampleSheet$Cell.type | is.na(pcaCellClassify)),],paste0(qcOutFolder, "SamplesPredictedDiffCellTypePCAOutlier.csv"))


predLabelledCellType<-sampleSheet$Cell.type == pcaCellClassify
predLabelledCellType[is.na(predLabelledCellType)]<-FALSE
sampleSheet$predLabelledCellType<-predLabelledCellType

## for time being take more stringent approach.
add.gdsn(gfile, 'QCdata', val = sampleSheet, replace = TRUE)

## need to close gds file in order to open in another R session
closefn.gds(gfile)

save(betas.scores, mahDistPCA, file = paste0(qcOutFolder,"WithinCellPCAValues.rdata"))

## add outcome to qc summary

sampleSheet<-sampleSheet[rownames(QCSum),]
QCSum<-cbind(QCSum, sampleSheet$predLabelledCellType)
colnames(QCSum)[ncol(QCSum)]<-"predLabelledCellType"

write.csv(QCSum, paste0(qcOutFolder,"PassQCStatusAllSamples.csv"))