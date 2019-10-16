## two normalisation strategies 
## 1) across all samples
## 2) within cell type
## script also performs sample filtering ## to add

library(bigmelon)
setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE, allow.fork = TRUE)

normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)

load(qcData)
load("QCmetrics/CellTypePCA.rdata")

## create new gfile with only samples that pass QC
if(exists(normgdsFile)){
	file.remove(normgdsFile) ## delete if already exists
}
normfile <- createfn.gds(normgdsFile)
for(node in ls.gdsn(gfile)){
	copyto.gdsn(node = normfile, source = index.gdsn(gfile, node), name = node)
}

## filter samples
passQC<-QCmetrics$intensPASS & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex)
write.csv(QCmetrics[!passQC,], "QCmetrics/ExcludedSamples.csv")
add.gdsn(gfile, 'QCoutcome', val = passQC, replace = TRUE)

## exclude cell type outliers
passCellQC<-rep("FALSE", nrow(QCmetrics))
passCellQC[which(passQC == TRUE)]<-!as.logical(cellOutlier)

add.gdsn(gfile, 'CellQCoutcome', val = passCellQC, replace = TRUE)
subSet(normfile, i=1:length(rownames(normfile)), j=(1:length(colnames(normfile)))[passCellQC]) ## i is probes, j is samples
## replace pData with QC metrics
QCmetrics<-QCmetrics[passCellQC,]
add.gdsn(normfile, 'pData', val = QCmetrics, replace = TRUE)


## repeat detection p value filtering
pfilter(normfile)

## need to close gds file in order to open in another R session
closefn.gds(gfile)

## normalise all samples together
dasen(normfile, node="normbeta")
normbetas<-index.gdsn(normfile, "normbeta")[,]

## need to extract below to run normalisation on each cell type
meth<-methylated(normfile)[,]
unmeth<-unmethylated(normfile)[,]
probeType<-fData(normfile)$Type
rawbetas<-betas(normfile)[,]

cellTypes<-unique(pData(normfile)$Cell.type)
projects<-unique(pData(normfile)$Project)

celltypeNormbeta<-matrix(NA, nrow = nrow(meth), ncol = ncol(meth))
rownames(celltypeNormbeta)<-rownames(rawbetas)
colnames(celltypeNormbeta)<-colnames(rawbetas)
for(entry in projects){
	for(each in cellTypes){
		index<-which(pData(normfile)$Cell.type == each & pData(normfile)$Project == entry)
		if(length(index) > 5){
			## perform normalisation within cell type & study
			## compare effect of normalisation within cell type
			celltypeNormbeta[,index]<-dasen(meth[,index], unmeth[,index], probeType)
		}
	}
}

add.gdsn(normfile, 'celltypenormbeta', val = celltypeNormbeta, replace = TRUE)

closefn.gds(normfile)