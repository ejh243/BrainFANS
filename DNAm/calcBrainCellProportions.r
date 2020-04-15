source("FunctionsForBrainCellProportionsPrediction.r")

library(bigmelon)
library(pheatmap)
library(glmnet)

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
QCSum<-read.csv(paste0(qcOutFolder,"PassQCStatusAllSamples.csv"), stringsAsFactors = FALSE, row.names = 1)

passQC<-rownames(QCSum)[which(QCSum$predLabelledCellType == "TRUE")]
QCmetrics<-read.gdsn(index.gdsn(gfile, "QCdata"))
QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

rawbetas<-betas(gfile)[,]
rawbetas<-rawbetas[,match(passQC, colnames(rawbetas))]

cellTypes<-unique(QCmetrics$Cell.type)
cellTypes<-cellTypes[!is.na(cellTypes)]
## sort so colours lines up correctly
cellTypes<-sort(cellTypes)
col_pal<-c("darkgreen", "darkblue", "darkmagenta", "deeppink", "darkgray") ## assumes celltypes are order alphabetically


load("RefDataForCellCompEstimation.rdata")
counts.all <- projectCellType(rawbetas[rownames(braincelldata), ], braincelldata)

## plot results
pdf(paste0(qcOutFolder, "PredictedCellComposition.pdf"),width = 12)
par(mfrow = c(2,3))
for(each in cellTypes){
	barplot(t(counts.all[which(QCmetrics$Cell.type == each),])*100, col = col_pal, main = each, ylab = "% estimated", names.arg = rep("", sum(QCmetrics$Cell.type == each)))
}
plot(0,1, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("center", colnames(counts.all), col = col_pal[1:ncol(counts.all)], pch = 15)

par(mfrow = c(2,2))
for(i in 1:ncol(counts.all)){
	boxplot(counts.all[,i]*100 ~ QCmetrics$Cell.type, col = col_pal, ylab = paste("%", colnames(counts.all)[i], "estimated"), xlab = "Cell type")
}

#aggregate(counts , by = list(QCmetrics$Cell.type), mean)

## filter out outlier samples

counts.sub<-counts.all[which(QCmetrics$predLabelledCellType == "TRUE"),]
sample.sub<-QCmetrics[which(QCmetrics$predLabelledCellType == "TRUE"),]

## plot results
par(mfrow = c(2,3))
for(each in unique(sample.sub$Cell.type)){
	barplot(t(counts.sub[which(sample.sub$Cell.type == each),])*100, col = col_pal, main = each, ylab = "% estimated")
}

par(mfrow = c(2,2))
for(i in 1:ncol(counts.all)){
	boxplot(counts.sub[,i]*100 ~ sample.sub$Cell.type, col = col_pal, ylab = paste("%", colnames(counts.all)[i], "estimated"), xlab = "Cell type")
}

#aggregate(counts ~ QCmetrics$Cell.type, mean)
dev.off()


#################

### repeat cellullar composition estmation excluding double negative

#################


load("RefDataForCellCompEstimationNoDoubleNeg.rdata")
counts.nodneg <- projectCellType(rawbetas[rownames(braincelldata), ], braincelldata)

## plot results
pdf(paste0(qcOutFolder, "PredictedCellCompositionNoDoubleNeg.pdf"),width = 12)
par(mfrow = c(2,3))
for(each in cellTypes){
	barplot(t(counts.nodneg[which(QCmetrics$Cell.type == each),])*100, col = col_pal[-1], main = each, ylab = "% estimated", names.arg = rep("", sum(QCmetrics$Cell.type == each)))
}
plot(0,1, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("center", colnames(counts.nodneg), col = col_pal[c(1:ncol(counts.nodneg))+1], pch = 15)

par(mfrow = c(2,2))
for(i in 1:ncol(counts.nodneg)){
	boxplot(counts.nodneg[,i]*100 ~ QCmetrics$Cell.type, col = col_pal, ylab = paste("%", colnames(counts.nodneg)[i], "estimated"), xlab = "Cell type")
}

#aggregate(counts , by = list(QCmetrics$Cell.type), mean)

## filter out outlier samples

counts.sub<-counts.nodneg[which(QCmetrics$predLabelledCellType == "TRUE"),]
sample.sub<-QCmetrics[which(QCmetrics$predLabelledCellType == "TRUE"),]

## plot results
par(mfrow = c(2,3))
for(each in unique(sample.sub$Cell.type)){
	barplot(t(counts.sub[which(sample.sub$Cell.type == each),])*100, col = col_pal[-1], main = each, ylab = "% estimated",names.arg = rep("", sum(sample.sub$Cell.type == each)))
}
plot(0,1, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("center", colnames(counts.nodneg), col = col_pal[c(1:ncol(counts.nodneg))+1], pch = 15)

par(mfrow = c(2,2))
for(i in 1:ncol(counts.nodneg)){
	boxplot(counts.sub[,i]*100 ~ sample.sub$Cell.type, col = col_pal, ylab = paste("%", colnames(counts.nodneg)[i], "estimated"), xlab = "Cell type")
}

#aggregate(counts ~ QCmetrics$Cell.type, mean)
dev.off()

## compare output with and without double negative
pdf(paste0(qcOutFolder, "ComparePredictedCellCompositionWithoutDoubleNeg.pdf"),width = 12, height = 4)
par(mfrow = c(1,3))
for(each in c("IRF8","NeuN","Sox10")){
	plot(counts.all[,each], counts.nodneg[,each], pch = 16,main = each, xlab = "With DNeg", ylab = "Without DNeg", col = col_pal[as.factor(QCmetrics$Cell.type)])
	abline(a = 0, b = 1)
}
dev.off()

