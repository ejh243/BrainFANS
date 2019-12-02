
## use QC metrics for phenotype data!!
## within adult samples do PCA to cluster sample types
thresBS<-80

library(bigmelon)
library(pheatmap)
library(glmnet)

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
load(qcData)
passQC<-QCmetrics$Basename[QCmetrics$intensPASS & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex)]

sampleSheet<-QCmetrics
sampleSheet<-sampleSheet[match(passQC, sampleSheet$Basename),]

rawbetas<-gfile[,, node = "betas"]
rawbetas<-rawbetas[,match(passQC, colnames(rawbetas))]
auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
rawbetas<-rawbetas[auto.probes,]

sampleKeep<-which(sampleSheet$Project == "MRC")
rawbetas<-rawbetas[,sampleKeep]
sampleSheet<-sampleSheet[sampleKeep,]
sampleSheet$Cell.type<-as.factor(as.character(sampleSheet$Cell.type))


cellTypes<-unique(sampleSheet$Cell.type)
cellCols<-rainbow(length(cellTypes))[as.factor(sampleSheet$Cell.type)]

## filter out NAs
rawbetas<-na.omit(rawbetas)

## filter out SNPs
rawbetas<-rawbetas[-grep("rs", rownames(rawbetas)),]

sample_anno<-sampleSheet[,c("Age","Sex", "Cell.type")]
sample_anno$Age<-as.numeric(sample_anno$Age)
rownames(sample_anno)<-sampleSheet$Basename
sigma<-apply(rawbetas, 1, sd)

## initial cluster
pheatmap(rawbetas[order(sigma, decreasing = TRUE)[1:500],], annotation_col = sample_anno,  show_colnames = FALSE, show_rownames = FALSE, cutree_cols = 4, main = "Pre Sample Type Filtering")


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

mahDistPCA<-matrix(data = NA, ncol = ncol(cellMedPCA), nrow = nrow(sampleSheet))
colnames(mahDistPCA)<-colnames(cellMedPCA)
for(each in colnames(cellMedPCA)){
	mahDistPCA[,each]<-mahalanobis(betas.scores, cellMedPCA[,each], cov_sigmaPCA[[each]], na.rm = TRUE)
}

y_lim<-range(mahDistPCA)
par(mfrow = c(2,2))
for(each in colnames(cellMedPCA)){
	boxplot(mahDistPCA[,each] ~ sampleSheet$Cell.type, main = paste("Comparision with ", each), col = rainbow(4), ylab = "Mahalanobis distance", xlab = "Labelled cell type")
}
closestCellTypePCA<-colnames(mahDistPCA)[unlist(apply(mahDistPCA, 1, which.min))]

table(sampleSheet$Cell.type,closestCellTypePCA)

###

pheatmap(rawbetas[order(sigma, decreasing = TRUE)[1:500],which(sampleSheet$Cell.type == closestCellTypePCA)], annotation_col = sample_anno[which(sampleSheet$Cell.type == closestCellTypePCA),],  show_colnames = FALSE, show_rownames = FALSE, cutree_cols = 4, main = "Mahalanobis Distance of PCs")

## look at PCA plot


predLabelledCellType<-sampleSheet$Cell.type == closestCellTypePCA

sampleSheet<-cbind(sampleSheet, closestCellTypePCA, predLabelledCellType)
write.csv(sampleSheet[which(predLabelledCellType == "FALSE"),], "QCmetrics/SamplesPredictedDiffCellTypePCAMahDist.csv")


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
for(i in 1:4){
	pcaCellClassify[which(betas.scores[,1] < upperBound[1,i] & betas.scores[,1] > lowerBound[1,i] & betas.scores[,2] < upperBound[2,i] & betas.scores[,2] > lowerBound[2,i])]<-colnames(lowerBound)[i]
}

sampleSheet<-cbind(sampleSheet, pcaCellClassify)

write.csv(sampleSheet[which(pcaCellClassify !=  sampleSheet$Cell.type | is.na(pcaCellClassify)),],"QCmetrics/SamplesPredictedDiffCellTypePCAOutlier.csv")

pheatmap(rawbetas[order(sigma, decreasing = TRUE)[1:500],which(sampleSheet$Cell.type == pcaCellClassify)], annotation_col = sample_anno[which(sampleSheet$Cell.type == pcaCellClassify),],  show_colnames = FALSE, show_rownames = FALSE, cutree_cols = 4, main = "PCA Outliers")

### plot of PCA to compare methods

par(mfrow = c(1,3))
plot(betas.scores[,1], betas.scores[,2], xlab = "PC 1", ylab = "PC 2", col = as.factor(sampleSheet$Cell.type))
legend("topright", pch = 16, col = palette()[1:4], colnames(lowerBound))
plot(betas.scores[,1], betas.scores[,2], pch = c(4,1)[as.factor(sampleSheet$Cell.type == closestCellTypePCA)], xlab = "PC 1", ylab = "PC 2", col = as.factor(sampleSheet$Cell.type), main = "Method 1")
points(cellMedPCA[1,], cellMedPCA[2,], pch = 16, col = as.factor(colnames(cellMedPCA)))

plot(betas.scores[,1], betas.scores[,2], xlab = "PC 1", ylab = "PC 2", col = as.factor(sampleSheet$Cell.type), pch = c(4,1)[as.factor(sampleSheet$Cell.type == pcaCellClassify)], main = "Method 2")
## add in NAs
points(betas.scores[is.na(pcaCellClassify),1], betas.scores[is.na(pcaCellClassify),2],col = as.factor(sampleSheet$Cell.type[is.na(pcaCellClassify)]), pch = 3)

for(i in 1:4){
	polygon(c(lowerBound[1,i], lowerBound[1,i], upperBound[1,i], upperBound[1,i]), c(lowerBound[2,i],upperBound[2,i],upperBound[2,i],lowerBound[2,i] ), border = palette()[i])
}

## compare filtering of both methods:

table(sampleSheet$Cell.type == closestCellTypePCA, sampleSheet$Cell.type == pcaCellClassify & !is.na(pcaCellClassify))

predLabelledCellType<-sampleSheet$Cell.type == pcaCellClassify
predLabelledCellType[is.na(predLabelledCellType)]<-FALSE
sampleSheet$predLabelledCellType<-predLabelledCellType
## for time being take more stringent approach.

add.gdsn(gfile, 'QCdata', val = sampleSheet, replace = TRUE)

##save all QC data

## need to close gds file in order to open in another R session
closefn.gds(gfile)