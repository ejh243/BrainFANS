##---------------------------------------------------------------------#
##
## Title: Comparision of normalisation strategies
##
## Purpose of script: compare effects of normalisation across all samples and within cell type.
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# uses wateRmelon metrics to compare normalisation strategies

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
set.seed(1234)
args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refPath <- args[2]

gdsFile <-file.path(dataDir, "2_gds/raw.gds")
normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)
qcOutFolder<-file.path(dataDir, "/2_gds/QCmetrics")
resPath<-file.path(dataDir, "4_analysis/methodsDevelopment/")

cellTypes <- c("Double-", "NeuN+", "Sox10+")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(bigmelon)
library(RColorBrewer)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

cellCols<-brewer.pal(n = length(cellTypes), name = "Paired")

setwd(dataDir)


gfile<-openfn.gds(normgdsFile, readonly = FALSE)
rawbetas<-betas(gfile)[,]

# need to extract below to run normalisation
meth<-methylated(gfile)[,]
unmeth<-unmethylated(gfile)[,]
celltypeNormbeta<-index.gdsn(gfile, "celltypenormbeta")[,]
#dasen(gfile, node="normbeta")
#normbeta<-index.gdsn(gfile, "normbeta")[,]
closefn.gds(gfile)

# list of samples that pass QC
QCSum<-read.csv(file.path(dataDir, "/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv"), row.names = 1, stringsAsFactors = FALSE)
passQC<-QCSum$Basename[which(QCSum$passQCS3)]
rm(QCSum)

QCmetrics<-read.csv(paste0(qcOutFolder,"/QCMetricsPostCellTypeClustering.csv"), stringsAsFactors = FALSE)
QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

## filter to NeuN+, Double- and Sox10 only
QCmetrics <-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]
rawbetas<-rawbetas[,QCmetrics$Basename]
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
meth<-meth[,QCmetrics$Basename]
unmeth<-unmeth[,QCmetrics$Basename]


probeAnnot<-read.table(file.path(refPath, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
probeAnnot<-probeAnnot[match(rownames(rawbetas), probeAnnot$probeID),]

# nprmalise all together
normbeta<-dasen(meth, unmeth, probeAnnot$designType)
rm(meth, unmeth)

## remove NAs 
rawbetas<-na.omit(rawbetas)
normbeta<-na.omit(normbeta)
celltypeNormbeta<-na.omit(celltypeNormbeta)

## filter to common probe set
probesKeep<-intersect(rownames(rawbetas), rownames(normbeta))
probesKeep<-intersect(probesKeep, rownames(celltypeNormbeta))

rawbetas<-rawbetas[probesKeep,]
normbeta<-normbeta[probesKeep,]
celltypeNormbeta<-celltypeNormbeta[probesKeep,]

#----------------------------------------------------------------------#
# CALCULATE NORMALISATION METRICS
#----------------------------------------------------------------------#


matOut<-matrix(data = NA, nrow = 3,ncol = 3)
rownames(matOut)<-c("raw", "normTog", "normSep")
colnames(matOut)<-c("iDMR", "genki", "seabi")

if(length( grep("rs", rownames(rawbetas))) > 0){
	matOut[1,2]<-mean(genki(rawbetas))
}
if(length( grep("rs", rownames(normbeta))) > 0){
	matOut[2,2]<-mean(genki(normbeta))
}
if(length( grep("rs", rownames(celltypeNormbeta))) > 0){
	matOut[3,2]<-mean(genki(celltypeNormbeta))
}

idmr<-intersect(iDMR(), rownames(rawbetas))

matOut[1,1]<-dmrse_row(rawbetas, idmr)
matOut[2,1]<-dmrse_row(normbeta, idmr)
matOut[3,1]<-dmrse_row(celltypeNormbeta, idmr)
rm(idmr)

probeAnnot<-probeAnnot[match(probesKeep, probeAnnot$probeID),]
x.probes<-probeAnnot$chrm == "chrX"

# To speed it up sample a random subset
indexSites<-sample(1:nrow(normbeta), 200000)
matOut[1,3]<-seabi(rawbetas[indexSites,], sex = QCmetrics$Sex, X = x.probes[indexSites])
matOut[2,3]<-seabi(normbeta[indexSites,], sex = QCmetrics$Sex, X = x.probes[indexSites])
matOut[3,3]<-seabi(celltypeNormbeta[indexSites,], sex = QCmetrics$Sex, X = x.probes[indexSites])

write.csv(matOut, file.path(resPath, "CompareNormalisationStrategies.csv"))


#----------------------------------------------------------------------#
# CALCULATE SAMPLE LEVEL NORMALISATION METRIC
#----------------------------------------------------------------------#

qualDat.tog<-qual(rawbetas, normbeta)
qualDat.sep<-qual(rawbetas, celltypeNormbeta)


## look at normalisation effects
plotLim<-range(c(qualDat.tog$rmsd, qualDat.sep$rmsd))
pdf(file.path(resPath, "CompareNormalisationStrategies.pdf"), height = 6, width = 10)
par(mfrow = c(1,2))
boxplot(qualDat.tog$rmsd ~ QCmetrics$Cell.type, 
ylab = "root mean square error", main = "Normalised together", xlab = "Cell type", col = cellCols,
ylim = plotLim)
boxplot(qualDat.sep$rmsd ~ QCmetrics$Cell.type, 
ylab = "root mean square error", main = "Normalised separately", xlab = "Cell type", col = cellCols,
ylim = plotLim)
dev.off()


