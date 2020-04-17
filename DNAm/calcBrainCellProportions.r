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

#############
### First ref panel all Exeter FANS fractions including double negative
#############

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
dev.off()

#aggregate(counts , by = list(QCmetrics$Cell.type), mean)


#################
### Second ref panel all Exeter FANS fractions excluding double negative
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

dev.off()

## compare output with and without double negative
pdf(paste0(qcOutFolder, "ComparePredictedCellCompositionWithoutDoubleNeg.pdf"),width = 12, height = 4)
par(mfrow = c(1,3))
for(each in c("IRF8","NeuN","Sox10")){
	plot(counts.all[,each], counts.nodneg[,each], pch = 16,main = each, xlab = "With DNeg", ylab = "Without DNeg", col = col_pal[as.factor(QCmetrics$Cell.type)])
	abline(a = 0, b = 1)
}
dev.off()

##############
### Third EpiGABA FANS samples
##############

load("RefDataEpiGABAForCellCompEstimation.rdata")
counts.epigaba <- projectCellType(rawbetas[intersect(rownames(braincelldata), rownames(rawbetas)), ], braincelldata[intersect(rownames(braincelldata), rownames(rawbetas)),])

## plot results
pdf(paste0(qcOutFolder, "PredictedCellCompositionEpiGabaRef.pdf"),width = 12)
par(mfrow = c(2,3))
for(each in cellTypes){
	barplot(t(counts.epigaba[which(QCmetrics$Cell.type == each),])*100, col = c("red", "blue", "green"), main = each, ylab = "% estimated", names.arg = rep("", sum(QCmetrics$Cell.type == each)))
}
plot(0,1, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("center", colnames(counts.epigaba), col = c("red", "blue", "green"), pch = 15)

par(mfrow = c(2,2))
for(i in 1:ncol(counts.epigaba)){
	boxplot(counts.epigaba[,i]*100 ~ QCmetrics$Cell.type, col = c("red", "blue", "green"), ylab = paste("%", colnames(counts.epigaba)[i], "estimated"), xlab = "Cell type")
}

dev.off()

pdf(paste0(qcOutFolder, "CompareCellCompositionEpiGabaBDRRef.pdf"),width = 10, height = 13)
par(mfrow = c(ncol(counts.all), ncol(counts.epigaba)))
for(i in 1:ncol(counts.all)){
	for(j in 1:ncol(counts.epigaba)){
		plot(counts.all[,i], counts.epigaba[,j], pch = 16, xlab = colnames(counts.all)[i], ylab = colnames(counts.epigaba)[j], col = col_pal[as.factor(QCmetrics$Cell.type)])
		abline(a = 0, b = 1)
	}
}
dev.off()

##############
### Fourth ref panel Merge EpiGABA and Exeter FANS samples
##############

load("RefDataForCellCompEstimationAll.rdata")
counts.merge <- projectCellType(rawbetas[intersect(rownames(braincelldata), rownames(rawbetas)), ], braincelldata[intersect(rownames(braincelldata), rownames(rawbetas)),])

## plot results
pdf(paste0(qcOutFolder, "PredictedCellCompositionMergedRef.pdf"),width = 12)
par(mfrow = c(2,3))
for(each in cellTypes){
	barplot(t(counts.merge[which(QCmetrics$Cell.type == each),])*100, col = c("red", "blue", "green", col_pal[2:4]), main = each, ylab = "% estimated", names.arg = rep("", sum(QCmetrics$Cell.type == each)))
}
plot(0,1, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("center", colnames(counts.merge), col = c("red", "blue", "green", col_pal[2:4]), pch = 15)

par(mfrow = c(2,3))
for(i in 1:ncol(counts.merge)){
	boxplot(counts.merge[,i]*100 ~ QCmetrics$Cell.type, col = c("red", "blue", "green", col_pal[2:4]), ylab = paste("%", colnames(counts.merge)[i], "estimated"), xlab = "Cell type")
}

dev.off()

##############
### Fifth ref c("GABAneurons","GLIA","GLUneurons","Sox10","IRF8")
### Fourth ref panel Merge EpiGABA and Exeter FANS samples
##############

load("RefDataForCellCompEstimationAllNoNeuN.rdata")
counts.merge4 <- projectCellType(rawbetas[intersect(rownames(braincelldata), rownames(rawbetas)), ], braincelldata[intersect(rownames(braincelldata), rownames(rawbetas)),])

## plot results
pdf(paste0(qcOutFolder, "PredictedCellCompositionMergedRefNoNeuN.pdf"),width = 12)
par(mfrow = c(2,3))
for(each in cellTypes){
	barplot(t(counts.merge4[which(QCmetrics$Cell.type == each),])*100, col = c("red", "blue", "green", col_pal[c(2,4)]), main = each, ylab = "% estimated", names.arg = rep("", sum(QCmetrics$Cell.type == each)))
}
plot(0,1, type = "n", xlab = "", ylab = "", axes = FALSE)
legend("center", colnames(counts.merge4), col = c("red", "blue", "green", col_pal[c(2,4)]), pch = 15)

par(mfrow = c(2,3))
for(i in 1:ncol(counts.merge4)){
	boxplot(counts.merge4[,i]*100 ~ QCmetrics$Cell.type, col = c("red", "blue", "green", col_pal[c(2,4)]), ylab = paste("%", colnames(counts.merge4)[i], "estimated"), xlab = "Cell type")
}

dev.off()

## compare with and without NeuN+ve
pdf(paste0(qcOutFolder, "ComparePredictedCellCompositionWithoutNeuN.pdf"),width = 12, height = 8)
par(mfrow = c(2,3))
for(each in c("GABAneurons","GLIA","GLUneurons","Sox10","IRF8")){
	plot(counts.merge[,each], counts.merge4[,each], pch = 16,main = each, xlab = "With NeuN", ylab = "Without NeuN", col = col_pal[as.factor(QCmetrics$Cell.type)])
	abline(a = 0, b = 1)
}
dev.off()