## two normalisation strategies 
## 1) across all samples
## 2) within cell type


args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

gdsFile <-file.path(dataDir, "/2_gds/raw.gds")
qcOutFolder<-file.path(dataDir, "/2_gds/QCmetrics")
normData<-file.path(dataDir, "/3_normalised/normalised.rdata")

library(bigmelon)

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
QCSum<-read.csv(file.path(dataDir, "/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv"), row.names = 1, stringsAsFactors = FALSE)
passQC<-QCSum$Basename[which(QCSum$passQCS3)]


QCmetrics<-read.csv(paste0(qcOutFolder,"/QCMetricsPostCellTypeClustering.csv"), stringsAsFactors = FALSE)
QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)

# create new gfile with only samples that pass QC
if(exists(normgdsFile)){
	file.remove(normgdsFile) ## delete if already exists
}
normfile <- createfn.gds(normgdsFile)
for(node in ls.gdsn(gfile)){
	copyto.gdsn(node = normfile, source = index.gdsn(gfile, node), name = node)
}

## filter out samples that fail QC
## need to match to basename not colnames!!!
rawbetas<-betas(normfile)[,]
subSet(normfile, i=1:length(rownames(normfile)), j=match(passQC, colnames(rawbetas))) 

## need to close gds file in order to open in another R session
closefn.gds(gfile)


## repeat detection p value filtering
#pfilter(normfile)

## normalise all samples together
dasen(normfile, node="normbeta")

## need to extract below to run normalisation on each cell type
meth<-methylated(normfile)[,]
unmeth<-unmethylated(normfile)[,]
probeType<-fData(normfile)$Type
rawbetas<-betas(normfile)[,]


cellTypes<-unique(QCmetrics$Cell.type)

celltypeNormbeta<-matrix(NA, nrow = nrow(meth), ncol = ncol(meth))
rownames(celltypeNormbeta)<-rownames(rawbetas)
colnames(celltypeNormbeta)<-colnames(rawbetas)
for(each in cellTypes){
	index<-which(QCmetrics$Cell.type == each)
	if(length(index) > 2){
		## perform normalisation within cell type & study
		## compare effect of normalisation within cell type
		celltypeNormbeta[,index]<-dasen(meth[,index], unmeth[,index], probeType)
	}
}


add.gdsn(normfile, 'celltypenormbeta', val = celltypeNormbeta, replace = TRUE)

closefn.gds(normfile)

save(celltypeNormbeta, QCmetrics, file = normData)