## after removing poor quality samples look within each cell type for additional outliers. 
## uses a PCA approach, only if there are at least 5 samples
## do this separately for two studies
## add some intital exclusion criteria

### parameters for sample filtering
thresBS<-80

library(bigmelon)
library(pheatmap)
library(glmnet)

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
load(qcData)
passQC<-QCmetrics$Basename[QCmetrics$intensPASS & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex)]

QCmetrics<-QCmetrics[QCmetrics$intensPASS & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex),]

studyInd<-QCmetrics$Project
cellTypes<-unique(QCmetrics$Cell.type)
cellCols<-rainbow(length(cellTypes))[as.factor(QCmetrics$Cell.type)]

rawbetas<-betas(gfile)[,]
rawbetas<-rawbetas[,match(QCmetrics$Basename, colnames(rawbetas))]
auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
rawbetas<-rawbetas[auto.probes,]

## first do heatmap to check how things look intially
sample_anno<-QCmetrics[,c("Cell.type", "Age","Sex")]
sample_anno$Age<-as.numeric(sample_anno$Age)
rownames(sample_anno)<-QCmetrics$Basename
sigma<-list()

cellPCA<-list()

cellOutlier<-rep("FALSE", length = ncol(rawbetas))

for(study in unique(studyInd)){
	for(each in cellTypes){
		sampleInd<-which(QCmetrics$Cell.type == each & QCmetrics$Project == study)
		if(length(sampleInd) > 5){ 
			## PCA within each cell type
			cellPCA[[study]][[each]]<-prcomp(t(na.omit(rawbetas[,sampleInd])))
			
			## denote anything more than 3 SD from mean an outlier on 1st PC
			mu<-mean(cellPCA[[study]][[each]]$x[,1], na.rm = TRUE)
			sigma<-sd(cellPCA[[study]][[each]]$x[,1], na.rm = TRUE)
			outlierindex<-which(cellPCA[[study]][[each]]$x[,1] > mu +3*sigma | cellPCA[[study]][[each]]$x[,1] < mu -3*sigma)
			if(length(outlierindex) > 0){
				cellOutlier[sampleInd[outlierindex]]<-"TRUE"
			}
		}
	}
}

save(cellPCA, cellOutlier, file = "QCmetrics/CellTypePCA.rdata")
	

## need to close gds file in order to open in another R session
closefn.gds(gfile)


