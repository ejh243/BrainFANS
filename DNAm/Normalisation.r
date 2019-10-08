## two normalisation strategies 
## 1) across all samples
## 2) within cell type
## script also performs sample filtering ## to add

library(bigmelon)
setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE, allow.fork = TRUE)

load(qcData)

## create new gfile with only samples that pass QC

if(exists('gdsFiles/norm.gds')){
	file.remove('gdsFiles/norm.gds') ## delete if already exists
}
normfile <- createfn.gds('gdsFiles/norm.gds')
for(node in ls.gdsn(gfile)){
	copyto.gdsn(node = normfile, source = index.gdsn(gfile, node), name = node)
}


## filter samples
passQC<-QCmetrics$M.median > intensThres & QCmetrics$U.median > intensThres & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) != as.character(QCmetrics$Sex)

write.csv(QCmetrics[!passQC,], "QCmetrics/ExcludedSamples.csv")

add.gdsn(gfile, 'QCoutcome', val = passQC, replace = TRUE)
subSet(normfile, i=1:length(rownames(normfile)), j=(1:length(colnames(normfile)))[passQC]) # If outliers is a vector of logical values where TRUE = an outlier ## i is probes, j is samples

## detection p value filtering
pfilter(normfile)

## need to close gds file in order to open in another R session
closefn.gds(gfile)

## normalise all samples together
dasen(normfile, node="normbeta")
normbetas<-index.gdsn(normfile, "normbeta")[,]

## filter to adult samples only 


## need to extract below to run normalisation on each cell type
meth<-methylated(normfile)[,]
unmeth<-unmethylated(normfile)[,]
probeType<-fData(normfile)$Type
rawbetas<-betas(normfile)[,]

cellTypes<-unique(pData(normfile)$Cell.type)

celltypeNormbeta<-matrix(NA, nrow = nrow(meth), ncol = ncol(meth))
rownames(celltypeNormbeta)<-rownames(rawbetas)
colnames(celltypeNormbeta)<-colnames(rawbetas)
for(each in cellTypes){
    if(length(which(pData(normfile)$Cell.type == each)) > 5){
		## perform normalisation within cell type
		## compare effect of normalisation within cell type
		celltypeNormbeta[,which(pData(normfile)$Cell.type == each)]<-dasen(meth[,which(pData(normfile)$Cell.type == each)], unmeth[,which(pData(normfile)$Cell.type == each)], probeType)
	}
}

add.gdsn(normfile, 'celltypenormbeta', val = celltypeNormbeta, replace = TRUE)

closefn.gds(normfile)