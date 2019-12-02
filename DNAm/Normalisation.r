## two normalisation strategies 
## 1) across all samples
## 2) within cell type
## script also performs sample filtering ## to add
thresBS<-80
source("") ## enter config file
library(bigmelon)
setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
#load(qcData)
QCmetrics<-read.gdsn(index.gdsn(gfile, "QCdata"))
passQC<-QCmetrics$Basename[as.logical(QCmetrics$intensPASS) & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex) & as.logical(QCmetrics$predLabelledCellType)]
write.csv(QCmetrics[-match(passQC, QCmetrics$Basename),], "QCmetrics/ExcludedSamples.csv")

QCmetrics<-QCmetrics[as.logical(QCmetrics$intensPASS) & QCmetrics$bisulfCon > thresBS & as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex) & as.logical(QCmetrics$predLabelledCellType),]



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
pfilter(normfile)

## normalise all samples together
dasen(normfile, node="normbeta")

## need to extract below to run normalisation on each cell type
meth<-methylated(normfile)[,]
unmeth<-unmethylated(normfile)[,]
probeType<-fData(normfile)$Type
rawbetas<-betas(normfile)[,]


cellTypes<-unique(QCmetrics$Cell.type)
projects<-unique(QCmetrics$Project)

celltypeNormbeta<-matrix(NA, nrow = nrow(meth), ncol = ncol(meth))
rownames(celltypeNormbeta)<-rownames(rawbetas)
colnames(celltypeNormbeta)<-colnames(rawbetas)
for(entry in projects){
	for(each in cellTypes){
		index<-which(QCmetrics$Cell.type == each & QCmetrics$Project == entry)
		if(length(index) > 5){
			## perform normalisation within cell type & study
			## compare effect of normalisation within cell type
			celltypeNormbeta[,index]<-dasen(meth[,index], unmeth[,index], probeType)
		}
	}
}

add.gdsn(normfile, 'celltypenormbeta', val = celltypeNormbeta, replace = TRUE)

closefn.gds(normfile)