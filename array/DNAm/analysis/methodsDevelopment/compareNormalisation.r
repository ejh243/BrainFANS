##---------------------------------------------------------------------#
##
## Title: Comparision of normalisation strategies
##
## Purpose of script: compare effects of normalisation across all samples and within cell type.
##
## Author: Eilis Hannon
##
## Date Created: 2020
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# uses wateRmelon metrics to compre normalised datasets

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
gdsFile <-file.path(dataDir, "2_gds/raw.gds")
normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)
qcOutFolder<-file.path(dataDir, "/2_gds/QCmetrics")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(bigmelon)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

setwd(dataDir)

gfile<-openfn.gds(normgdsFile, readonly = FALSE)

normbetas<-index.gdsn(gfile, "normbeta")[,]
celltypeNormbeta<-index.gdsn(gfile, "celltypenormbeta")[,]
rawbetas<-betas(gfile)[,]

## remove NAs 
rawbetas<-na.omit(rawbetas)
normbetas<-na.omit(normbetas)
celltypeNormbeta<-na.omit(celltypeNormbeta)

QCmetrics<-read.csv(paste0(qcOutFolder,"/QCMetricsPostCellTypeClustering.csv"), stringsAsFactors = FALSE)
QCmetrics<-QCmetrics[match(colnames(rawbetas), QCmetrics$Basename),]

#----------------------------------------------------------------------#
# CALCULATE NORMALISATION METRICS
#----------------------------------------------------------------------#


matOut<-matrix(data = NA, nrow = 3,ncol = 3)
rownames(matOut)<-c("raw", "normTog", "normSep")
colnames(matOut)<-c("iDMR", "genki", "seabi")

if(length( grep("rs", rownames(rawbetas))) > 0){
	matOut[1,2]<-mean(genki(rawbetas))
}
if(length( grep("rs", rownames(normbetas))) > 0){
	matOut[2,2]<-mean(genki(normbetas))
}
if(length( grep("rs", rownames(celltypeNormbeta))) > 0){
	matOut[3,2]<-mean(genki(celltypeNormbeta))
}


## filter to common probes
probes<-intersect(intersect(rownames(rawbetas), rownames(normbetas)), rownames(celltypeNormbeta))
rawbetas<-rawbetas[probes,]
normbetas<-normbetas[probes,]
celltypeNormbeta<-celltypeNormbeta[probes,]

probeAnno<-fData(gfile)
probeAnno<-probeAnno[probes,]
x.probes<-probeAnno$chr == "chrX"


idmr<-intersect(iDMR(), rownames(rawbetas))

matOut[1,1]<-dmrse_row(rawbetas, idmr)
matOut[2,1]<-dmrse_row(normbetas, idmr)
matOut[3,1]<-dmrse_row(celltypeNormbeta, idmr)

matOut[1,3]<-seabi(rawbetas, sex = QCmetrics$Sex, X = x.probes)
matOut[2,3]<-seabi(normbetas, sex = QCmetrics$Sex, X = x.probes)
matOut[3,3]<-seabi(celltypeNormbeta, sex = QCmetrics$Sex, X = x.probes)

write.csv(matOut, paste0(qcOutFolder, "CompareNormalisationStrategies.csv"))

pheno<-QCmetrics

#----------------------------------------------------------------------#
# CALCULATE SAMPLE LEVEL NORMALISATION METRIC
#----------------------------------------------------------------------#

qualDat.tog<-qual(rawbetas, normbetas)
qualDat.sep<-qual(rawbetas, celltypeNormbeta)

cellCols<-c("orange", "darkblue", "darkmagenta", "deeppink", "darkgray") ## assumes celltypes are order alphabetically


## look at normalisation effects
pdf(paste0(qcOutFolder, "CompareNormalisationStrategies.pdf"), height = 6, width = 10)
par(mfrow = c(1,2))
boxplot(qualDat.tog$rmsd ~ QCmetrics$Cell.type, ylab = "root mean square error", main = "Normalised together", xlab = "Cell type")
boxplot(qualDat.sep$rmsd ~ QCmetrics$Cell.type, ylab = "root mean square error", main = "Normalised separately", xlab = "Cell type")

## look at distribution of beta values
densityPlot(normbetas, sampGroups = QCmetrics$Cell.type,pal = cellCols)
densityPlot(celltypeNormbeta, sampGroups = QCmetrics$Cell.type)
dev.off()


## reextract cell norm corrected betas to save in bespoke DNAm data object
celltypeNormbeta<-index.gdsn(gfile, "celltypeNormbeta")[,]
closefn.gds(gfile)


## remove cross-hybridising & snp probes
crosshyb<-read.table(paste(refFiles, "/EPICArray/CrossHydridisingProbes_McCartney.txt", sep = ""), stringsAsFactors = FALSE)
snpProbes<-read.table(paste(refFiles, "/EPICArray/SNPProbes_McCartney.txt", sep = ""), stringsAsFactors = FALSE, header = TRUE)
crosshyb2<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM1.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes2<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM4.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes3<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM5.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes4<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM6.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes<-snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01),]
snpProbes2<-snpProbes2[which(snpProbes2$AF > 0.01),]
snpProbes3<-snpProbes3[which(snpProbes3$AF > 0.01),]
snpProbes4<-snpProbes4[which(snpProbes4$AF > 0.01),]

dist<-cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist<-apply(dist, 1, min)
snpProbes4<-snpProbes4[which(dist <=10),]

remove<-intersect(rownames(celltypeNormbeta),unique(c(crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE)))

celltypeNormbeta<-celltypeNormbeta[!rownames(celltypeNormbeta) %in% remove,]

save(celltypeNormbeta, pheno, file = normData)

