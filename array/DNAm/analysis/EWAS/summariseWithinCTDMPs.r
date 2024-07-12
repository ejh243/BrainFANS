## ---------------------------------------------------------------------#
##
## Title: Classify & Compare EWAS Results from regression within each cell type
##
## Purpose of script: Characterise EWAS results and produce summary plots
##
## ---------------------------------------------------------------------#

library(qqman)

args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refPath <- args[2]

normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resPath <- file.path(dataDir, "4_analysis/EWAS")

cellTypes <- c("Double-", "NeuN+", "Sox10+")

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
probeAnnot$chrm <- gsub("chr", "", probeAnnot$chrm)
probeAnnot$start <- probeAnnot$start+1

for(i in 1:3){
    res[[i]]<- cbind(res[[i]], probeAnnot[, c("chrm", "start", "GeneNames", "GeneClasses", "CGI", "CGIPosition")])
    res[[i]]<- res[[i]][which(res[[i]][,"chrm"] != "Y"),]
}

#----------------------------------------------------------------------#
# SUMMARY PLOTS
#----------------------------------------------------------------------#

## QQ plots

tiff(file.path(resPath, "Plots","QQplotsWithinCTLRSCZEffectNoCCNoCETYGO.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
for(i in 1:3){
    qq(res[[i]][,"NullModel_SCZ_P"], main = cellTypes[i])
}
dev.off()

tiff(file.path(resPath, "Plots","QQplotsWithinCTLRSCZEffectNoCCWithCETYGO.tiff"), width = 1200, height = 400)
par(mfrow = c(1,3))
for(i in 1:3){
    qq(res[[i]][,"FullModel_SCZ_P"], main = cellTypes[i])
}
dev.off()

## Compare with and without CETYGO & cell composition

par(mfrow = c(3,2))
for(i in 1:3){
    if("CCModel_SCZ_P" %in% colnames(res[[i]])){
        plot(-log10(res[[i]][,"FullModel_SCZ_P"]), -log10(res[[i]][,"CCModel_SCZ_P"]), pch = 16, 
        xlab = "With CETYGO", ylab = "Without CETYGO", main = cellTypes[i])
        abline(a = 0, b = 1)
        plot(-log10(res[[i]][,"CCModel_SCZ_P"]), -log10(res[[i]][,"NullModel_SCZ_P"]), pch = 16, 
        xlab = "With cell composition", ylab = "Without cell composition", main = cellTypes[i])
        abline(a = 0, b = 1)
    } else {
        plot(-log10(res[[i]][,"FullModel_SCZ_P"]), -log10(res[[i]][,"NullModel_SCZ_P"]), pch = 16, 
        xlab = "With CETYGO", ylab = "Without CETYGO", main = cellTypes[i])
        abline(a = 0, b = 1)
    }
}

# Manhattan plots

jpeg(file.path(resPath, "Plots", "ManhattanPlotsWithinCTLM.jpeg"), width = 600, height = 600, quality = 100)
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
for(i in 1:3){
    res[[i]][,"chrm"][which(res[[i]][,"chrm"] == "X")]<-"23"
    res[[i]][,"chrm"]<-as.numeric(as.character(res[[i]][,"chrm"]))

    resMan<-data.frame("SNP" = rownames(res[[i]]), "P" = res[[i]][,"FullModel_SCZ_P"], "CHR" = res[[i]][,"chrm"],"BP" = res[[i]][,"start"])
    resMan<-na.omit(resMan)
    manhattan(resMan, genomewideline = -log10(9e-8), suggestiveline = -log10(5e-5), main = cellTypes[i])
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
cbind()

write.csv(nDMPs, file.path(resPath, "DMPCountsPerCelltypeAcrossWithinCTModels.csv"))