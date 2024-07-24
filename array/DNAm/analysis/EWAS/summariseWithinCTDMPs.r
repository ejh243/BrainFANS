## ---------------------------------------------------------------------#
##
## Title: Classify & Compare EWAS Results from regression within each cell type
##
## Purpose of script: Characterise EWAS results and produce summary plots
##
## ---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refPath <- args[2]

normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resPath <- file.path(dataDir, "4_analysis/EWAS")

cellTypes <- c("Double-", "NeuN+", "Sox10+")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(qqman)
library(tidyr)
library(ggplot2)
library(ggpubr)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

res<-list()
for (each in cellTypes) {
    load(file.path(resPath, paste0(each, "LM.rdata")))
    res[[each]] <- outtab
}
rm(outtab)

#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#

## filter out cross hyb & SNP probes
crosshyb <- read.table(file.path(refPath, "CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter <- read.csv(file.path(refPath, "EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes <- read.table(file.path(refPath, "SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
crosshyb2 <- read.csv(file.path(refPath, "Pidsley_SM1.csv"), stringsAsFactors = FALSE)
snpProbes2 <- read.csv(file.path(refPath, "Pidsley_SM4.csv"), stringsAsFactors = FALSE)
snpProbes3 <- read.csv(file.path(refPath, "Pidsley_SM5.csv"), stringsAsFactors = FALSE)
snpProbes4 <- read.csv(file.path(refPath, "Pidsley_SM6.csv"), stringsAsFactors = FALSE)

snpProbes <- snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01), ]
snpProbes2 <- snpProbes2[which(snpProbes2$AF > 0.01), ]
snpProbes3 <- snpProbes3[which(snpProbes3$AF > 0.01), ]
snpProbes4 <- snpProbes4[which(snpProbes4$AF > 0.01), ]

dist <- cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist <- apply(dist, 1, min)
snpProbes4 <- snpProbes4[which(dist <= 10), ]

remove <- unique(c(tofilter$IlmnID, crosshyb[, 1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))
remove <- intersect(remove, rownames(res[[1]]))

if (length(remove) > 0) {
    for (i in 1:3) {
        res[[i]] <- res[[i]][-match(remove, rownames(res[[i]])), ]
    }
}


# add gene annotation
probeAnnot <- read.table(file.path(refPath, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
probeAnnot <- probeAnnot[match(rownames(res[[1]]), probeAnnot$probeID), ]
probeAnnot[["chrm"]] <- gsub("chr", "", probeAnnot[["chrm"]], fixed = TRUE)
probeAnnot$start <- probeAnnot$start+1

for(i in 1:3){
    res[[i]]<- cbind(res[[i]], probeAnnot[, c("chrm", "start", "GeneNames", "GeneClasses", "CGI", "CGIPosition")])
    res[[i]]<- res[[i]][which(res[[i]][,"chrm"] != "Y"),]
}

#----------------------------------------------------------------------#
# SUMMARY PLOTS
#----------------------------------------------------------------------#

## QQ plots

png(file.path(resPath, "Plots","QQplotsWithinCTLRSCZEffect.png"), width = 1200, height = 1200, res = 200)
par(mfrow = c(3,3))
par(mar = c(4,4,1,1))
par(mgp = c(2,0.75,0))
for(i in 1:3){
    qq(res[[i]][,"NullModel_SCZ_P"], main = cellTypes[i])
}
for(i in 1:3){
    if("CCModel_SCZ_P" %in% colnames(res[[i]])){
        qq(res[[i]][,"CCModel_SCZ_P"], main = cellTypes[i])
    } else {
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    }
}
for(i in 1:3){
    qq(res[[i]][,"FullModel_SCZ_P"], main = cellTypes[i])
}

dev.off()

## Compare with and without CETYGO & cell composition
png(file.path(resPath, "Plots","ScatterplotsWithinCTLRSCZEffect.png"), width = 800, height = 1200, res = 200)
par(mfrow = c(3,2))
par(mar = c(4,4,1,1))
par(mgp = c(2,0.75,0))
for(i in 1:3){
    if("CCModel_SCZ_P" %in% colnames(res[[i]])){
        plot(-log10(res[[i]][,"FullModel_SCZ_P"]), -log10(res[[i]][,"CCModel_SCZ_P"]), pch = 16, 
        xlab = "With CETYGO", ylab = "Without CETYGO", main = cellTypes[i], col = "darkgrey")
        abline(a = 0, b = 1)
        plot(-log10(res[[i]][,"CCModel_SCZ_P"]), -log10(res[[i]][,"NullModel_SCZ_P"]), pch = 16, 
        xlab = "With cell composition", ylab = "Without cell composition", main = cellTypes[i], col = "darkgrey")
        abline(a = 0, b = 1)
    } else {
        plot(-log10(res[[i]][,"FullModel_SCZ_P"]), -log10(res[[i]][,"NullModel_SCZ_P"]), pch = 16, 
        xlab = "With CETYGO", ylab = "Without CETYGO", main = cellTypes[i], col = "darkgrey")
        abline(a = 0, b = 1)
    }
}
dev.off()

# Manhattan plots

png(file.path(resPath, "Plots", "ManhattanPlotsWithinCTLM.png"), width = 1400, height = 1000, res = 200)
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
par(mgp = c(2,0.75,0))
for(i in 1:3){
    res[[i]][,"chrm"][which(res[[i]][,"chrm"] == "X")]<-"23"
    res[[i]][,"chrm"]<-as.numeric(as.character(res[[i]][,"chrm"]))

    resMan<-data.frame("SNP" = rownames(res[[i]]), "P" = res[[i]][,"FullModel_SCZ_P"], "CHR" = res[[i]][,"chrm"],"BP" = res[[i]][,"start"])
    resMan<-na.omit(resMan)
    manhattan(resMan, genomewideline = -log10(9e-8), suggestiveline = -log10(1e-5), main = cellTypes[i])
} 
dev.off()


#----------------------------------------------------------------------#
# IDENTIFY DMPS
#----------------------------------------------------------------------#

## count DMPs at diff thresholds
nDMPs<-matrix(ncol = 9, nrow = 4)
colnames(nDMPs)<-outer(c("FullModel", "CCModel", "NullModel"),  cellTypes, paste, sep = ":")
rowIndex<-1
for(thres in c(9e-8, 1e-7, 1e-6, 1e-5)){
    for(i in 1:3){
        nDMPs[rowIndex,(i*3)-2]<-length(which(res[[i]][,"FullModel_SCZ_P"] < thres))
        if("CCModel_SCZ_P" %in% colnames(res[[i]])){
            nDMPs[rowIndex,(i*3)-1]<-length(which(res[[i]][,"CCModel_SCZ_P"] < thres))
        }
        nDMPs[rowIndex,(i*3)]<-length(which(res[[i]][,"NullModel_SCZ_P"] < thres))
    }
    rowIndex<-rowIndex+1
}

write.csv(nDMPs, file.path(resPath, "Tables", "DMPCountsPerCelltypeAcrossWithinCTModels.csv"))

for(each in cellTypes){
	write.csv(res[[each]][which(res[[each]][,"FullModel_SCZ_P"] < thres),], file.path(resPath, "Tables", paste0("DiscoveryDMPsWithin", each, "LM.csv")))
}

## COMPARE DMPS ACROSS CELL TYPES
thres<-1e-5

dmpList<-NULL
for(i in 1:3){
    dmpList<-rbind(dmpList, cbind(rownames(res[[i]])[which(res[[i]][,"FullModel_SCZ_P"] < thres)], cellTypes[i]))
}
colnames(dmpList)<-c("ProbeID", "DiscoveryCellType")



dmpRes<-cbind(dmpList, res[[1]][dmpList[,1], c("FullModel_SCZ_P", "FullModel_SCZ_coeff", "FullModel_SCZ_SE")],
res[[2]][dmpList[,1], c("FullModel_SCZ_P", "FullModel_SCZ_coeff", "FullModel_SCZ_SE")],
res[[3]][dmpList[,1], c("FullModel_SCZ_P", "FullModel_SCZ_coeff", "FullModel_SCZ_SE")])
colnames(dmpRes)<-c("ProbeID", "DiscoveryCellType", outer(c("P", "Coeff", "SE"), cellTypes, paste, sep = ":"))



diffLong<-pivot_longer(data.frame(dmpRes[,c(2,4,7,10)]), cols = gsub("\\+|-", "\\.", paste("Coeff", cellTypes, sep = ":")))
diffLong[["name"]] <- gsub("Coeff\\.", "", diffLong[["name"]], fixed = TRUE)
diffLong$value<-abs(diffLong$value)
pdf(file.path(resPath, "Plots","ViolinPlotCelltypeEffectsDiscoveryDMPsLMWithinCTs.pdf"), width = 8, height = 4)
ggplot(diffLong, aes(x = name, y = value, fill = name)) + geom_violin() +
xlab("Cell Type") + ylab("Mean difference") + 
stat_summary(fun=mean, geom="point", size=2, color="black") + facet_wrap (~DiscoveryCellType)
dev.off()

plotLim<-range(dmpRes[,c(4,7,10)])

plotByStatus <- function(dmpCellType, constrastCellType, dmpRes, plotLim){
    subset(dmpRes, DiscoveryCellType == dmpCellType) %>%
    ggplot(aes(x = `Coeff:Double-`, y = `Coeff:NeuN+`)) + geom_point() +
    xlab(paste0(dmpCellType, "ve")) + ylab(paste0(dmpCellType, "ve")) + ggtitle(paste0(dmpCellType, "ve DMPs")) +
    xlim(plotLim) + ylim(plotLim) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0)

}

p1<- plotByStatus("NeuN+", "Double-", dmpRes, plotLim)
p2<- plotByStatus("NeuN+", "Sox10+", dmpRes, plotLim)
p3<- plotByStatus("Sox10+", "Double-", dmpRes, plotLim)
p4<- plotByStatus("Sox10+", "NeuN+", dmpRes, plotLim)
p5<- plotByStatus("Double-", "NeuN+", dmpRes, plotLim)
p6<- plotByStatus("Double-", "Sox10+", dmpRes, plotLim)


pdf(file.path(resPath, "Plots","ScatterplotsCelltypeEffectsDiscoveryDMPsLMWithinCTs.pdf"), width = 12, height = 8)
ggarrange(p1, p3,p5, p2,p4,   p6, nrow = 2, ncol = 3)
dev.off()

load(file.path(resPath, "MLM.rdata"))
	
outtab<-cbind(outtab, outtab[,"SCZ_coeff"],
outtab[,"SCZ_coeff"]+outtab[,"NeuN_SCZ_coeff"],
outtab[,"SCZ_coeff"]+outtab[,"SOX10_SCZ_coeff"])
colnames(outtab)[(ncol(outtab)-2):ncol(outtab)]<-c("DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")


dmpRes<-cbind(dmpRes, outtab[dmpList[,"ProbeID"],
c("SCZ_P","CellType_SCZ_P","NeuN_SCZ_P","SOX10_SCZ_P","CellType_P", "DNeg_Mean_Diff", "NEUN_Mean_Diff", "SOX10_Mean_Diff")])

ggplot(dmpRes, aes(x = DiscoveryCellType, y = -log10(CellType_SCZ_P), fill = DiscoveryCellType)) + geom_violin() + 
stat_summary(fun=mean, geom="point", size=2, color="black") +
xlab("Discovery Cell Type") + ylab("Heterogeneity -log10P") + geom_hline(yintercept = -log10(0.05)) 


# To rule out power effects plot coeff against p-value

topSites<-NULL
for(i in 1:3){
    index<-order(res[[i]][,"FullModel_SCZ_P"])[1:100]
    topSites<-rbind(topSites, cbind(cellTypes[i], res[[i]][index,c("FullModel_SCZ_P" , "FullModel_SCZ_coeff")]))
}
colnames(topSites)[1]<-"cellType"

topSites[,2]<- -log10(topSites[,2])
topSites[,3] <- abs(topSites[,3])

ggplot(topSites, aes(x = FullModel_SCZ_P, y = FullModel_SCZ_coeff, color = cellType)) + geom_point() +
xlab("-log10(P)") + ylab("Mean difference")


## Plot DMPs
pos <- position_dodge(0.9)

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

dmpRes <- cbind(dmpRes, probeAnnot[match(rownames(dmpRes), probeAnnot$probeID), c("chrm", "start", "GeneNames", "GeneClasses", "CGI", "CGIPosition")])


for(each in cellTypes){
	subRes <- dmpRes[which(dmpRes[,"DiscoveryCellType"] == each),]
	pdf(file.path(resPath, "Plots", paste0("ViolinPlotDiscoveryDMPsWithin", each, "Models.pdf")), width = 6, height = 5)

	for(i in 1:nrow(subRes)){
		tmpDat<-data.frame("CellType" = QCmetrics$Cell.type, "Phenotype" = QCmetrics$Phenotype, "DNAm" = celltypeNormbeta[rownames(subRes)[i],])
		p<-ggplot(tmpDat, aes(x=CellType, y=DNAm, fill = Phenotype)) + 
	  geom_violin(position = pos, scale = 'width') + ggtitle(paste(rownames(subRes)[i], subRes$GeneNames[i])) +
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos)
		print(p)
	}
	dev.off()
}
