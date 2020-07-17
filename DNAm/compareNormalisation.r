
## compare effects of normalisation across all samples and within cell type.
## use metrics from original wateRmelon paper to do comparision

args<-commandArgs(trailingOnly = TRUE)
source(args[1])

library(bigmelon)

setwd(dataDir)
normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)
gfile<-openfn.gds(normgdsFile, readonly = FALSE)

normbetas<-index.gdsn(gfile, "normbeta")[,]
celltypenormbeta<-index.gdsn(gfile, "celltypenormbeta")[,]
rawbetas<-betas(gfile)[,]

## remove NAs 
rawbetas<-na.omit(rawbetas)
normbetas<-na.omit(normbetas)
celltypenormbeta<-na.omit(celltypenormbeta)

QCmetrics<-read.gdsn(index.gdsn(gfile, "QCdata"))
QCmetrics<-QCmetrics[match(colnames(rawbetas), QCmetrics$Basename),]

## apply wateRmelon metrics to compare normalisation
matOut<-matrix(data = NA, nrow = 3,ncol = 3)
rownames(matOut)<-c("raw", "normTog", "normSep")
colnames(matOut)<-c("iDMR", "genki", "seabi")

if(length( grep("rs", rownames(rawbetas))) > 0){
	matOut[1,2]<-mean(genki(rawbetas))
}
if(length( grep("rs", rownames(normbetas))) > 0){
	matOut[2,2]<-mean(genki(normbetas))
}
if(length( grep("rs", rownames(celltypenormbeta))) > 0){
	matOut[3,2]<-mean(genki(celltypenormbeta))
}


## filter to common probes
probes<-intersect(intersect(rownames(rawbetas), rownames(normbetas)), rownames(celltypenormbeta))
rawbetas<-rawbetas[probes,]
normbetas<-normbetas[probes,]
celltypenormbeta<-celltypenormbeta[probes,]

probeAnno<-fData(gfile)
probeAnno<-probeAnno[probes,]
x.probes<-probeAnno$chr == "chrX"


idmr<-intersect(iDMR(), rownames(rawbetas))

matOut[1,1]<-dmrse_row(rawbetas, idmr)
matOut[2,1]<-dmrse_row(normbetas, idmr)
matOut[3,1]<-dmrse_row(celltypenormbeta, idmr)

matOut[1,3]<-seabi(rawbetas, sex = QCmetrics$Sex, X = x.probes)
matOut[2,3]<-seabi(normbetas, sex = QCmetrics$Sex, X = x.probes)
matOut[3,3]<-seabi(celltypenormbeta, sex = QCmetrics$Sex, X = x.probes)

write.csv(matOut, paste0(qcOutFolder, "CompareNormalisationStrategies.csv"))

pheno<-QCmetrics


qualDat.tog<-qual(rawbetas, normbetas)
qualDat.sep<-qual(rawbetas, celltypenormbeta)

cellCols<-c("orange", "darkblue", "darkmagenta", "deeppink", "darkgray") ## assumes celltypes are order alphabetically


## look at normalisation effects
pdf(paste0(qcOutFolder, "CompareNormalisationStrategies.pdf"), height = 6, width = 10)
par(mfrow = c(1,2))
boxplot(qualDat.tog$rmsd ~ QCmetrics$Cell.type, ylab = "root mean square error", main = "Normalised together", xlab = "Cell type")
boxplot(qualDat.sep$rmsd ~ QCmetrics$Cell.type, ylab = "root mean square error", main = "Normalised separately", xlab = "Cell type")

## look at distribution of beta values
densityPlot(normbetas, sampGroups = QCmetrics$Cell.type,pal = cellCols)
densityPlot(celltypenormbeta, sampGroups = QCmetrics$Cell.type)
dev.off()


## reextract cell norm corrected betas to save in bespoke DNAm data object
celltypenormbeta<-index.gdsn(gfile, "celltypenormbeta")[,]
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

remove<-intersect(rownames(celltypenormbeta),unique(c(crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE)))

celltypenormbeta<-celltypenormbeta[!rownames(celltypenormbeta) %in% remove,]

save(celltypenormbeta, pheno, file = normData)

