##---------------------------------------------------------------------#
##
## Title: Determine power in other cell types to confirm cell specificity
##
## Purpose of script: Calculate power for typical DMP effect sizes in 
## other cell types to confirm cell specific enrichment of findings 
##
##---------------------------------------------------------------------#

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

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pwr)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

res<-list()
for (each in cellTypes) {
    load(file.path(resPath, paste0(each, "LM.rdata")))
    res[[each]] <- outtab
}
rm(outtab)


setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

sampleNumbers<-table(QCmetrics$Phenotype, QCmetrics$Cell.type)


#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#


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
    res[[i]]<- res[[i]][!res[[i]][,"chrm"] %in% c("Y", "*"),]
}



#----------------------------------------------------------------------#
# CALCULATE POWER
#----------------------------------------------------------------------#

## identify NeuN DMPs
pThres<-1e-6
DMPindex<-which(res[["NeuN+"]][,"NullModel_SCZ_P"] < pThres)
powerCalcs<-matrix(data = NA, nrow = length(DMPindex), ncol = 3)
colnames(powerCalcs)<-cellTypes
rownames(powerCalcs)<-rownames(res[["NeuN+"]])[DMPindex]
for(probe in rownames(res[["NeuN+"]])[DMPindex]){
    for(CT in cellTypes){
        sigma<-sd(celltypeNormbeta[probe,which(QCmetrics$Cell.type == CT)])
        powerCalcs[probe,CT] <- pwr.t2n.test(d = res[["NeuN+"]][probe,"NullModel_SCZ_coeff"]/sigma, 
            sig.level = 1e-6, 
            alternative = "two.sided", 
            n1 = sampleNumbers[1,CT], n2 = sampleNumbers[2,CT])$power
    }
}

powerCalcsLong<-pivot_longer(as.data.frame(powerCalcs), cols = everything())

pdf(file.path(resPath, "Plots","ViolinPlotPowerStatisticsDiscoveryNeuNDMPs.pdf"), width = 4, height = 4)
ggplot(powerCalcsLong, aes(x = name, y = value, fill = name)) + geom_violin() +
xlab("Cell Type") + ylab("Power") + 
stat_summary(fun=mean, geom="point", size=2, color="black") +
labs(fill = "Cell Type")
dev.off()



## BUT we would expect that cell specific effects for other cell types are at other sites
typicalES <- median(abs(res[["NeuN+"]][DMPindex,"NullModel_SCZ_coeff"]))
binSize <- 500
powerCalcs <- matrix(data = NA, nrow = binSize, ncol = length(cellTypes))
colnames(powerCalcs)<-cellTypes
for(each in cellTypes){
  # subset to cell type
  betasCell <- as.matrix(celltypeNormbeta[,QCmetrics$Cell.type == each])
  # calculate SD
  betasSD <- apply(as.matrix(betasCell),1,sd)
    betasSD <- as.data.frame(betasSD) %>% mutate(points_bin = ntile(betasSD, n=binSize))

    for(j in unique(betasSD$points_bin)){
      meanSD <- mean(betasSD[which(betasSD$points_bin == j),1])
      powerCalcs[j, each] <- pwr.t2n.test(d = typicalES/meanSD, 
            sig.level = 1e-6, 
            alternative = "two.sided", 
            n1 = sampleNumbers[1,CT], n2 = sampleNumbers[2,CT])$power
    }
}

pdf(file.path(resPath, "Plots","ViolinPlotPowerStatisticsAllSitesTypicalNeuNEffect.pdf"), width = 4, height = 4)
powerCalcsLong<-pivot_longer(as.data.frame(powerCalcs), cols = everything())
ggplot(powerCalcsLong, aes(x = name, y = value, fill = name)) + geom_violin() +
xlab("Cell Type") + ylab("Power") + 
stat_summary(fun=mean, geom="point", size=2, color="black")+
labs(fill = "Cell Type")
dev.off()

    
