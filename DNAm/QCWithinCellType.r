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
passQC<-QCmetrics$intensPASS & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex)

sampleSheet<-pData(gfile)[passQC,]
studyInd<-sampleSheet$Study
cellTypes<-unique(sampleSheet$Cell.type)
cellCols<-rainbow(length(cellTypes))[as.factor(sampleSheet$Cell.type)]

rawbetas<-betas(gfile)[,passQC]
auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
rawbetas<-rawbetas[auto.probes,]

## first do heatmap to check how things look intially
sample_anno<-sampleSheet[,c("Cell.type", "Age","Sex")]
sample_anno$Age<-as.numeric(sample_anno$Age)
rownames(sample_anno)<-sampleSheet$Basename
sigma<-list()

cellPCA<-list()

for(study in unique(studyInd)){
	for(each in cellTypes){
		sampleInd<-which(sampleSheet$Cell.type == each & sampleSheet$Study == study)
		if(length(sampleInd) > 5){ 
			## PCA within each cell type
			cellPCA[[study]][[each]]<-prcomp(t(na.omit(rawbetas[,sampleInd])))
		}
	}
}

save(cellPCA, file = "QCmetrics/CellTypePCA.rdata")
	

## need to close gds file in order to open in another R session
closefn.gds(gfile)
