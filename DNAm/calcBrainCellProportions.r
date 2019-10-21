source("FunctionsForBrainCellProportionsPrediction.r")

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

rawbetas<-betas(gfile)[,]
rawbetas<-rawbetas[,match(passQC, colnames(rawbetas))]
auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
rawbetas<-rawbetas[auto.probes,]

sampleKeep<-which(sampleSheet$Project == "MRC")
rawbetas<-rawbetas[,sampleKeep]
sampleSheet<-sampleSheet[sampleKeep,]
sampleSheet$Cell.type<-as.factor(as.character(sampleSheet$Cell.type))


load("RefDataForCellCompEstimation.rdata")
counts <- projectCellType(rawbetas[rownames(braincelldata), ], braincelldata)

## plot results
par(mfrow = c(2,2))
for(each in unique(sampleSheet$Cell.type)){
	barplot(t(counts[which(sampleSheet$Cell.type == each),])*100, col = brewer.pal(3, "Set1"), main = each, ylab = "% estimated", names.arg = sampleSheet$Samples.ID[which(sampleSheet$Cell.type == each)])
}

par(mfrow = c(1,3))
for(i in 1:ncol(counts)){
	boxplot(counts[,i]*100 ~ sampleSheet$Cell.type, col = brewer.pal(4, "Set1"), ylab = paste("%", colnames(counts)[i], "estimated"), xlab = "Cell type")
}

aggregate(counts , by = list(sampleSheet$Cell.type), mean)

## filter out outlier samples

counts.sub<-counts[which(sampleSheet$predLabelledCellType == "TRUE"),]
sample.sub<-sampleSheet[which(sampleSheet$predLabelledCellType == "TRUE"),]

## plot results
par(mfrow = c(2,2))
for(each in unique(sample.sub$Cell.type)){
	barplot(t(counts.sub[which(sample.sub$Cell.type == each),])*100, col = brewer.pal(3, "Set1"), main = each, ylab = "% estimated")
}

par(mfrow = c(1,3))
for(i in 1:ncol(counts)){
	boxplot(counts.sub[,i]*100 ~ sample.sub$Cell.type, col = brewer.pal(4, "Set1"), ylab = paste("%", colnames(counts)[i], "estimated"), xlab = "Cell type")
}

aggregate(counts ~ sampleSheet$Cell.type, mean)



