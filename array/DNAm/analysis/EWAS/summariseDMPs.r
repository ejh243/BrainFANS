##---------------------------------------------------------------------#
##
## Title: Classify & Compare EWAS Results
##
## Purpose of script: Characterise EWAS results and produce summary plots
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE CHARACTERISATION FUNCTIONS
#----------------------------------------------------------------------#

calcCTDiffs<-function(res){
	# use output of interaction model to calculate mean case control difference for each cell type 
	res<-cbind(res, res[,"SCZ_coeff"],
	res[,"SCZ_coeff"]+res[,"NeuN_SCZ_coeff"],
	res[,"SCZ_coeff"]+res[,"SOX10_SCZ_coeff"])
	colnames(res)[(ncol(res)-2):ncol(res)]<-c("DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")
	return(res)
}

calcMaxCTDiff<-function(res){
 res<-cbind(res, apply(abs(as.matrix(res[,c("DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")])), 1, max))
 colnames(res)[ncol(res)]<-"CTMaxDiff"
 return(res)
}


classifyDMPs<-function(res, pThres1 = 5e-5, pThres2 = 0.05){
	# identify DMPs for each cell type
	# effect in dneg if SCZ_P signif
	dneg_dmp <- res[,"SCZ_P"] < pThres1
	# effect in neun if SCZ_P signif and neuN_SCZ_P not signif or if SCZ_P not signif and NeuN_SCZ_P signif
	neun_dmp <- res[,"SCZ_P"] < pThres1 & res[,"NeuN_SCZ_P"] > pThres2 | res[,"SCZ_P"] > pThres2 & res[,"NeuN_SCZ_P"] < pThres1
	# effect in sox10 if SCZ_P signif and SOX10_SCZ_P not signif or if SCZ_P not signif and SOX10_SCZ_P signif
	sox10_dmp <- res[,"SCZ_P"] < pThres1 & res[,"SOX10_SCZ_P"] > pThres2 | res[,"SCZ_P"] > pThres2 & res[,"SOX10_SCZ_P"] < pThres1

	res<-cbind(res, dneg_dmp, neun_dmp, sox10_dmp)
	colnames(res)[(ncol(res)-2):ncol(res)]<-c("DNeg_DMP", "NEUN_DMP", "SOX10_DMP")
	return(res)
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refPath <- args[2]

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis/EWAS")

pos <- position_dodge(0.9)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
cellTypes<-unique(QCmetrics$Cell.type)

# load results
load(file.path(resPath, "LM.rdata"))
res.lm<-outtab
load(file.path(resPath, "MLM.rdata"))
res.mlm<-outtab
load(file.path(resPath, "RobustClusteredErrors.rdata"))
res.rce<-outtab

load("Analysis/CellType/ANOVABrainCellTypes.rdata")

#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#

## filter out cross hyb & SNP probes
crosshyb<-read.table(file.path(refPath, "CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter<-read.csv(file.path(refPath, "EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes<-read.table(file.path(refPath, "SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
crosshyb2<-read.csv(file.path(refPath, "Pidsley_SM1.csv"), stringsAsFactors = FALSE)
snpProbes2<-read.csv(file.path(refPath, "Pidsley_SM4.csv"), stringsAsFactors = FALSE)
snpProbes3<-read.csv(file.path(refPath, "Pidsley_SM5.csv"), stringsAsFactors = FALSE)
snpProbes4<-read.csv(file.path(refPath, "Pidsley_SM6.csv"), stringsAsFactors = FALSE)

snpProbes<-snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01),]
snpProbes2<-snpProbes2[which(snpProbes2$AF > 0.01),]
snpProbes3<-snpProbes3[which(snpProbes3$AF > 0.01),]
snpProbes4<-snpProbes4[which(snpProbes4$AF > 0.01),]

dist<-cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist<-apply(dist, 1, min)
snpProbes4<-snpProbes4[which(dist <=10),]

remove<-unique(c(tofilter$IlmnID, crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))
remove<-intersect(remove, rownames(res.lm))

if(length(remove) > 0){
	res.lm<-res.lm[-match(remove, rownames(res.lm)),]
	res.mlm<-res.mlm[-match(remove, rownames(res.mlm)),]
	res.rce<-res.rce[-match(remove, rownames(res.rce)),]
}

#----------------------------------------------------------------------#
# ANNOTATE RESULTS
#----------------------------------------------------------------------#

# add columns with cell type mean differences
res.lm<-calcCTDiffs(res.lm)
res.mlm<-calcCTDiffs(res.mlm)
res.rce<-calcCTDiffs(res.rce)

res.lm<-calcMaxCTDiff(res.lm)
res.mlm<-calcMaxCTDiff(res.mlm)
res.rce<-calcMaxCTDiff(res.rce)

# add columns to indicate cell type dmps
res.lm<-classifyDMPs(res.lm)
res.mlm<-classifyDMPs(res.mlm)
res.rce<-classifyDMPs(res.rce)

# add gene annotation
annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))
newfData <- newfData[rownames(res.lm), ]

res.lm<-cbind(res.lm, newfData[,c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Regulatory_Feature_Group", "Regulatory_Feature_Name", "Relation_to_Island")])
res.mlm<-cbind(res.mlm, newfData[,c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Regulatory_Feature_Group", "Regulatory_Feature_Name", "Relation_to_Island")])
res.rce<-cbind(res.rce, newfData[,c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Regulatory_Feature_Group", "Regulatory_Feature_Name", "Relation_to_Island")])

#----------------------------------------------------------------------#
# WRITE TABLES
#----------------------------------------------------------------------#

dmps.lm<-res.lm[rowSums(as.matrix(res.lm[,c("DNeg_DMP", "NEUN_DMP", "SOX10_DMP")])) > 0,]
dmps.mlm<-res.mlm[rowSums(as.matrix(res.mlm[,c("DNeg_DMP", "NEUN_DMP", "SOX10_DMP")])) > 0,]
dmps.rce<-res.rce[rowSums(as.matrix(res.rce[,c("DNeg_DMP", "NEUN_DMP", "SOX10_DMP")])) > 0,]

write.csv(dmps.lm, file.path(resPath, "Tables", "DiscoveryDMPsAllCellTypesLM.csv"))
write.csv(dmps.mlm, file.path(resPath, "Tables", "DiscoveryDMPsAllCellTypesMLM.csv"))
write.csv(dmps.rce, file.path(resPath, "Tables", "DiscoveryDMPsAllCellTypesCRR.csv"))

#----------------------------------------------------------------------#
# SUMMARY PLOTS
#----------------------------------------------------------------------#

tiff(file.path(resPath, "Plots","QQplotsMainSCZEffect.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
qq(res.lm[,"SCZ_P"])
qq(res.mlm[,"SCZ_P"])
qq(res.rce[,"SCZ_P"])
dev.off()

tiff(file.path(resPath, "Plots","QQplotsSCZCTIntEffect.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
qq(res.lm[,"CellType_SCZ_P"])
qq(res.mlm[,"CellType_SCZ_P"])
qq(res.rce[,"CellType_SCZ_P"])
dev.off()

tiff(file.path(resPath, "Plots","QQplotsSCZNEUNCTIntEffect.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
qq(res.lm[,"NeuN_SCZ_P"])
qq(res.mlm[,"NeuN_SCZ_P"])
qq(res.rce[,"NeuN_SCZ_P"])
dev.off()

tiff(file.path(resPath, "Plots","QQplotsSCZSOX10CTIntEffect.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
qq(res.lm[,"SOX10_SCZ_P"])
qq(res.mlm[,"SOX10_SCZ_P"])
qq(res.rce[,"SOX10_SCZ_P"])
dev.off()

### compare P-values across models
tiff(file.path(resPath, "Plots","Scatterlotlog10SCZMainEffectP.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
plot(-log10(res.lm[,"SCZ_P"]), -log10(res.mlm[,"SCZ_P"]), pch = 16, xlab = "OLS", ylab = "Mixed effects", cex = 1.1)
abline(a = 0, b = 1)
plot(-log10(res.lm[,"SCZ_P"]), -log10(res.rce[,"SCZ_P"]), pch = 16, xlab = "OLS", ylab = "Robust clustered", cex = 1.1)
abline(a = 0, b = 1)
plot(-log10(res.mlm[,"SCZ_P"]), -log10(res.rce[,"SCZ_P"]), pch = 16, xlab = "Mixed effects", ylab = "Robust clustered", cex = 1.1)
abline(a = 0, b = 1)
dev.off()

## manhattan plots
res.lm[,"chr"][which(res.lm[,"chr"] == "chrX")]<-"chr23"
res.lm[,"chr"][which(res.lm[,"chr"] == "chrY")]<-"chr24"
res.lm[,"chr"]<-gsub("chr", "", res.lm[,"chr"])
res.lm[,"chr"]<-as.numeric(as.character(res.lm[,"chr"]))

res<-data.frame("SNP" = rownames(res.lm), "P" = res.lm[,"SCZ_P"], "CHR" = res.lm[,"chr"],"BP" = res.lm[,"pos"])
res<-na.omit(res)
jpeg(file.path(resPath, "Plots", "ManhattanPlotsLM.jpeg"), width = 600, height = 600, quality = 100)
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))

res<-data.frame("SNP" = rownames(res.lm), "P" = res.lm[,"NeuN_SCZ_P"], "CHR" = res.lm[,"chr"],"BP" = res.lm[,"pos"])
res<-na.omit(res)
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))
res<-data.frame("SNP" = rownames(res.lm), "P" = res.lm[,"SOX10_SCZ_P"], "CHR" = res.lm[,"chr"],"BP" = res.lm[,"pos"])
res<-na.omit(res)
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))
dev.off()

res.mlm[,"chr"][which(res.mlm[,"chr"] == "chrX")]<-"chr23"
res.mlm[,"chr"][which(res.mlm[,"chr"] == "chrY")]<-"chr24"
res.mlm[,"chr"]<-gsub("chr", "", res.mlm[,"chr"])
res.mlm[,"chr"]<-as.numeric(as.character(res.mlm[,"chr"]))

res<-data.frame("SNP" = rownames(res.mlm), "P" = res.mlm[,"SCZ_P"], "CHR" = res.mlm[,"chr"],"BP" = res.mlm[,"pos"])
res<-na.omit(res)
jpeg(file.path(resPath, "Plots", "ManhattanPlotsMLM.jpeg"), width = 600, height = 600, quality = 100)
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))

res<-data.frame("SNP" = rownames(res.mlm), "P" = res.mlm[,"NeuN_SCZ_P"], "CHR" = res.mlm[,"chr"],"BP" = res.mlm[,"pos"])
res<-na.omit(res)
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))
res<-data.frame("SNP" = rownames(res.mlm), "P" = res.mlm[,"SOX10_SCZ_P"], "CHR" = res.mlm[,"chr"],"BP" = res.mlm[,"pos"])
res<-na.omit(res)
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))
dev.off()

res.rce[,"chr"][which(res.rce[,"chr"] == "chrX")]<-"chr23"
res.rce[,"chr"][which(res.rce[,"chr"] == "chrY")]<-"chr24"
res.rce[,"chr"]<-gsub("chr", "", res.rce[,"chr"])
res.rce[,"chr"]<-as.numeric(as.character(res.rce[,"chr"]))

res<-data.frame("SNP" = rownames(res.rce), "P" = res.rce[,"SCZ_P"], "CHR" = res.rce[,"chr"],"BP" = res.rce[,"pos"])
res<-na.omit(res)
jpeg(file.path(resPath, "Plots", "ManhattanPlotsCRR.jpeg"), width = 600, height = 600, quality = 100)
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))

res<-data.frame("SNP" = rownames(res.rce), "P" = res.rce[,"NeuN_SCZ_P"], "CHR" = res.rce[,"chr"],"BP" = res.rce[,"pos"])
res<-na.omit(res)
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))
res<-data.frame("SNP" = rownames(res.rce), "P" = res.rce[,"SOX10_SCZ_P"], "CHR" = res.rce[,"chr"],"BP" = res.rce[,"pos"])
res<-na.omit(res)
manhattan(res, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5))
dev.off()


## COMPARE DMPS ACROSS CELL TYPES
library(tidyr)
pdf(file.path(resPath, "Plots","ViolinPlotCelltypeEffectsDiscoveryDMPsLM.pdf"), width = 6, height = 6)
diffLong<-gather(data.frame(dmps.lm[,c("DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")]))
diffLong$key<-gsub("_Mean_Diff", "", diffLong$key)
diffLong$value<-abs(diffLong$value)
ggplot(diffLong, aes(x = key, y = value, fill = key)) + geom_violin() +xlab("Cell Type") + ylab("Mean difference") + ggtitle(paste("Across", nrow(dmps.lm), "Discovery DMPs")) + stat_summary(fun=mean, geom="point", size=2, color="black")
dev.off()

pdf(file.path(resPath, "Plots","ViolinPlotCelltypeEffectsDiscoveryDMPsMLM.pdf"), width = 6, height = 6)
diffLong<-gather(data.frame(dmps.mlm[,c("DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")]))
diffLong$key<-gsub("_Mean_Diff", "", diffLong$key)
diffLong$value<-abs(diffLong$value)
ggplot(diffLong, aes(x = key, y = value, fill = key)) + geom_violin() +xlab("Cell Type") + ylab("Mean difference") + ggtitle(paste("Across", nrow(dmps.mlm), "Discovery DMPs")) + stat_summary(fun=mean, geom="point", size=2, color="black")
dev.off()

pdf(file.path(resPath, "Plots","ViolinPlotCelltypeEffectsDiscoveryDMPsCRR.pdf"), width = 6, height = 6)
diffLong<-gather(data.frame(dmps.rce[,c("DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")]))
diffLong$key<-gsub("_Mean_Diff", "", diffLong$key)
diffLong$value<-abs(diffLong$value)
ggplot(diffLong, aes(x = key, y = value, fill = key)) + geom_violin() +xlab("Cell Type") + ylab("Mean difference") + ggtitle(paste("Across", nrow(dmps.rce), "Discovery DMPs")) + stat_summary(fun=mean, geom="point", size=2, color="black")
dev.off()
 
 
pdf(file.path(resPath, "Plots","CompareCelltypeEffectsDiscoveryDMPsLM.pdf"), width = 12, height = 4)
par(mfrow = c(1,3))
par(mar = c(4.5,4.5,1,1))
plot(dmps.lm[,"DNeg_Mean_Diff"], dmps.lm[,"NEUN_Mean_Diff"], pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff NeuN+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dmps.lm[,"DNeg_Mean_Diff"], dmps.lm[,"SOX10_Mean_Diff"], pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dmps.lm[,"NEUN_Mean_Diff"], dmps.lm[,"SOX10_Mean_Diff"], pch = 16, xlab = "Mean diff NeuN+ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)
dev.off()

pdf(file.path(resPath, "Plots","CompareCelltypeEffectsDiscoveryDMPsMLM.pdf"), width = 12, height = 4)
par(mfrow = c(1,3))
par(mar = c(4.5,4.5,1,1))
plot(dmps.mlm[,"DNeg_Mean_Diff"], dmps.mlm[,"NEUN_Mean_Diff"], pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff NeuN+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dmps.mlm[,"DNeg_Mean_Diff"], dmps.mlm[,"SOX10_Mean_Diff"], pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dmps.mlm[,"NEUN_Mean_Diff"], dmps.mlm[,"SOX10_Mean_Diff"], pch = 16, xlab = "Mean diff NeuN+ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)
dev.off()

pdf(file.path(resPath, "Plots","CompareCelltypeEffectsDiscoveryDMPsCRR.pdf"), width = 12, height = 4)
par(mfrow = c(1,3))
par(mar = c(4.5,4.5,1,1))
plot(dmps.rce[,"DNeg_Mean_Diff"], dmps.rce[,"NEUN_Mean_Diff"], pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff NeuN+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dmps.rce[,"DNeg_Mean_Diff"], dmps.rce[,"SOX10_Mean_Diff"], pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dmps.rce[,"NEUN_Mean_Diff"], dmps.rce[,"SOX10_Mean_Diff"], pch = 16, xlab = "Mean diff NeuN+ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)
dev.off()

## PLOT DMPS

library(ggplot2)
# Basic violin plot


pdf(file.path(resPath, "Plots", "ViolinPlotDiscoveryDMPsRCEModel.pdf"), width = 6, height = 5)
for(i in order(dmps.rce[,"CTMaxDiff"], decreasing = TRUE)){
	tmpDat<-data.frame("CellType" = QCmetrics$Cell.type, "Phenotype" = QCmetrics$Phenotype, "DNAm" = celltypeNormbeta[rownames(dmps.rce)[i],])
	p<-ggplot(tmpDat, aes(x=CellType, y=DNAm, fill = Phenotype)) + 
  geom_violin(position = pos, scale = 'width') + ggtitle(paste(rownames(dmps.rce)[i], dmps.rce$UCSC_RefGene_Name[i])) +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos)
	print(p)
}
dev.off()


pdf(file.path(resPath, "Plots", "ViolinPlotDiscoveryDMPsMLMModel.pdf"), width = 6, height = 5)
for(i in order(dmps.mlm[,"CTMaxDiff"], decreasing = TRUE)){
	if(dmps.mlm[i,"CTMaxDiff"] > 0.04){
		tmpDat<-data.frame("CellType" = QCmetrics$Cell.type, "Phenotype" = QCmetrics$Phenotype, "DNAm" = celltypeNormbeta[rownames(dmps.mlm)[i],])
		p<-ggplot(tmpDat, aes(x=CellType, y=DNAm, fill = Phenotype)) + 
	  geom_violin(position = pos, scale = 'width') + ggtitle(paste(rownames(dmps.mlm)[i], dmps.mlm$UCSC_RefGene_Name[i])) +
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos)
		print(p)
	}
}
dev.off()


# characterise DMPs

par(mfrow = c(1,3))
plot(-log10(dmps.lm[,"SCZ_P"]),-log10(dmps.lm[,"NeuN_SCZ_P"]), pch = 16, xlab = "-log10 Main Effect P", ylab = "log10 NeuN * Interaction P")
plot(-log10(dmps.lm[,"SCZ_P"]),-log10(dmps.lm[,"CellType_P"]), pch = 16, xlab = "-log10 Main Effect P", ylab = "log10 Cell type P")

hist(-log10(dmps.lm[which(res.lm[,"SCZ_P"] < 9e-8),"CellType_SCZ_P"]), xlab = "SCZ x Cell type P")
abline(v = -log10(0.05))
hist(-log10(res.lm[which(res.lm[,"SCZ_P"] < 9e-8),"CellType_P"]), xlab = "Cell type P")
abline(v = -log10(9e-8), col = "red")

par(mfrow = c(1,2))
plot(-log10(res.mlm[,16]),-log10(res.mlm[,"CellType_P"]), pch = 16, xlab = "-log10 Main Effect P", ylab = "log10 Cell type P")
plot(-log10(res.mlm[,16]),-log10(res.mlm[,"CellType_SCZ_P"]), pch = 16, xlab = "-log10 Main Effect P", ylab = "log10 Interaction P")
hist(-log10(res.mlm[which(res.mlm[,16] < 9e-8),"CellType_P"]), xlab = "Cell type -logP", col = "gray")
abline(v = -log10(9e-8), col = "red")
hist(-log10(res.mlm[which(res.mlm[,16] < 9e-8),"CellType_SCZ_P"]), xlab = "SCZ x Cell type -logP", col = "gray")
abline(v = -log10(0.05))

par(mfrow = c(1,2))
plot(-log10(res.rce[,"SCZ_P"]),-log10(res.rce[,"CellType_SCZ_P"]), pch = 16, xlab = "-log10 Main Effect P", ylab = "log10 Interaction P")
plot(-log10(res.rce[,"SCZ_P"]),-log10(res.rce[,"CellType_P"]), pch = 16, xlab = "-log10 Main Effect P", ylab = "log10 Cell type P")
hist(-log10(res.rce[order(res.rce[,"SCZ_P"])[1:25],"CellType_P"]), xlab = "Cell type P", col = "gray")
abline(v = -log10(9e-8), col = "red")
hist(-log10(res.rce[order(res.rce[,"SCZ_P"])[1:25],"CellType_SCZ_P"]), xlab = "SCZ x Cell type P", col = "gray")
abline(v = -log10(0.05))


par(mfrow = c(1,3))
par(mar = c(4.5,4.5,1,1))
plot(dneg, neun, pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff NeuN+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(dneg, sox, pch = 16, xlab = "Mean diff Double-ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)

plot(neun, sox, pch = 16, xlab = "Mean diff NeuN+ve", ylab = "Mean diff Sox10+ve", cex = 1.5, cex.axis = 1.3, cex.lab = 1.3)
abline(v = 0)
abline(h = 0)
abline(a = 0, b = 1)





## are any scz dmps sig diff between irf8 and dneg?

ctOut<-ctOut[rownames(res.mlm),]

par(mfrow = c(1,2))
table(ctOut[which(res.mlm[,16] < 9e-8),6] < 9e-8)
plot(ctOut[which(res.mlm[,16] < 9e-8),5], res.mlm[which(res.mlm[,16] < 9e-8),14]*-1, pch = 16
, col = as.factor(ctOut[which(res.mlm[,16] < 9e-8),6] < 9e-8), xlab = "IRF8vsDNeg", ylab = "SCZDiff")
abline(v = 0)
abline(h = 0)

table(ctOut[order(res.rce[,"SCZ_P"])[1:25],6] < 9e-8)
plot(ctOut[order(res.rce[,"SCZ_P"])[1:25],5], res.rce[order(res.rce[,"SCZ_P"])[1:25],14]*-1, pch = 16,
col = as.factor(ctOut[order(res.rce[,"SCZ_P"])[1:25],6] < 9e-8), xlab = "IRF8vsDNeg", ylab = "SCZDiff")
abline(v = 0)
abline(h = 0)

pdf("Analysis/Schizophrenia/CompareModels/BoxplotTop25DMPsRCEModel.pdf", width = 6, height = 5)
for(i in order(res.rce[,"SCZ_P"])[1:25]){
	boxplot(celltypenormbeta[rownames(res.rce)[i],] ~ pheno$Phenotype*pheno$Cell.type, names = rep(c("SCZ", "CON"), 3)
	, xlab = "", ylab = "DNAm", main = "", col = c("#6699CC", "#FF6666"))
	abline(v = seq(2.5, 5, 2), lty = 2)
	mtext(side = 3, line = 0, adj = 0, "Double Negative")
	mtext(side = 3, line = 0, adj = 0.5, "NeuN+ve")
	mtext(side = 3, line = 0, adj = 1, "Sox10+ve")
	mtext(side = 3, line = 2, paste(rownames(res.rce)[i], newfData$UCSC_RefGene_Name[i]))
	mtext(side = 3, line = 1, adj = 1, paste("P =", signif(res.rce[i,14],3), "Int P =",signif(res.rce[i,"CellType_SCZ_P"],3)))
}
dev.off()


pdf("Analysis/Schizophrenia/CompareModels/BoxplotTop25DMPsMLMModel.pdf", width = 6, height = 5)
for(i in order(res.mlm[,16])[1:25]){
	boxplot(celltypenormbeta[rownames(res.mlm)[i],] ~ pheno$Phenotype*pheno$Cell.type, names = rep(c("SCZ", "CON"), 3)
	, xlab = "", ylab = "DNAm", main = "", col = c("#6699CC", "#FF6666"))
	abline(v = seq(2.5, 5, 2), lty = 2)
	mtext(side = 3, line = 0, adj = 0, "Double Negative")
	mtext(side = 3, line = 0, adj = 0.5, "NeuN+ve")
	mtext(side = 3, line = 0, adj = 1, "Sox10+ve")
	mtext(side = 3, line = 2, paste(rownames(res.mlm)[i], newfData$UCSC_RefGene_Name[i]))
	mtext(side = 3, line = 1, adj = 1, paste("P =", signif(res.mlm[i,16],3), "Int P =",signif(res.mlm[i,"CellType_SCZ_P"],3)))
}
dev.off()


## replot excluding sites with significant difference between irf8 and double negative

pdf("Analysis/Schizophrenia/CompareModels/BoxplotTop25DMPsRCEModel.pdf", width = 12, height = 5)
for(i in order(res.rce[,"SCZ_P"])[1:25]){
	boxplot(celltypenormbeta[rownames(res.rce)[i],] ~ pheno$Phenotype*pheno$Cell.type, names = rep(c("SCZ", "CON"), 3)
	, xlab = "", ylab = "DNAm", main = "", col = c("#6699CC", "#FF6666"))
	abline(v = seq(2.5, 5, 2), lty = 2)
	mtext(side = 3, line = 0, adj = 0, "Double Negative")
	mtext(side = 3, line = 0, adj = 0.5, "NeuN+ve")
	mtext(side = 3, line = 0, adj = 1, "Sox10+ve")
	mtext(side = 3, line = 2, paste(rownames(res.rce)[i], newfData$UCSC_RefGene_Name[i]))
	mtext(side = 3, line = 1, adj = 1, paste("P =", signif(res.rce[i,16],3), "Int P =",signif(res.rce[i,"CellType_SCZ_P"],3)))
}
dev.off()


pdf("Analysis/Schizophrenia/CompareModels/BoxplotDMPsMLMNoDiffIRF8Model.pdf", width = 12, height = 5)
for(i in which(res.mlm[,16] < 9e-8 & ctOut[,6] > 0.05)){
	boxplot(celltypenormbeta[rownames(res.mlm)[i],] ~ pheno$Phenotype*pheno$Cell.type, names = rep(c("SCZ", "CON"), 3)
	, xlab = "", ylab = "DNAm", main = "", col = c("#6699CC", "#FF6666"))
	abline(v = seq(2.5, 5, 2), lty = 2)
	mtext(side = 3, line = 0, adj = 0, "Double Negative")
	mtext(side = 3, line = 0, adj = 0.5, "NeuN+ve")
	mtext(side = 3, line = 0, adj = 1, "Sox10+ve")
	mtext(side = 3, line = 2, paste(rownames(res.mlm)[i], newfData$UCSC_RefGene_Name[i]))
	mtext(side = 3, line = 1, adj = 1, paste("P =", signif(res.mlm[i,16],3), "Int P =",signif(res.mlm[i,"CellType_SCZ_P"],3)))
}
dev.off()