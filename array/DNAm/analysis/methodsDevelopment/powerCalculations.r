##---------------------------------------------------------------------#
##
## Title: Perform power calculations for cell specific EWAS
##
## Purpose of script: Quantify power for within cell type EWAS under different
## experimental settings
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refPath <- args[2]

resPath<-file.path(dataDir, "4_analysis/methodsDevelopment/")
normData<-file.path(dataDir, "/3_normalised/normalised.rdata")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(BrainPower)
library(pwr)
library(doParallel)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

load(normData)
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# subset data to autosomes
man <- fread(file.path(refPath, "MethylationEPIC_v-1-0_B5.csv"), skip=7, fill=TRUE, data.table=F)

auto.probes<-man$IlmnID[man$CHR != "X" & man$CHR != "Y" & man$CHR != ""]
celltypeNormbeta<-celltypeNormbeta[row.names(celltypeNormbeta) %in% auto.probes,]

cellTypes<-unique(QCmetrics$Cell.type)

allSDs <- matrix(data = NA, nrow = nrow(celltypeNormbeta), ncol = length(cellTypes))
rownames(allSDs)<-rownames(celltypeNormbeta)
colnames(allSDs)<-cellTypes
for(each in cellTypes){
  # subset to cell type
  betasCell <- as.matrix(celltypeNormbeta[,QCmetrics$Cell.type == each])
  # calculate SD
  betasSD <- apply(as.matrix(betasCell),1,sd)
  allSDs[,each] <- betasSD
}

cellCols<-brewer.pal(n = ncol(allSDs), name = "Paired")
names(cellCols)<-colnames(allSDs)

#----------------------------------------------------------------------#
# PLOTS OF SD DISTRIBUTION
#----------------------------------------------------------------------#

plotdf <- as.data.frame(allSDs)
plotdf$ID <- row.names(plotdf)
plotdf <- melt(plotdf)
colnames(plotdf)[2:3] <- c("CellType", "SD")


# density plot
pdf(file.path(resPath, "densityPlotSDs.pdf"), width = 5, height = 5)
ggplot(plotdf, aes(x=SD, colour=CellType))+
  scale_color_manual(values = cellCols) +
  geom_density(size = 2)+
  xlab("SD (DNA methylation)")+
  theme(legend.key = element_rect(fill = NA))
dev.off()

# violin plot
pdf(file.path(resPath, "violinPlotSDs.pdf"), width = 5, height = 5)
ggplot(plotdf, aes(x=CellType, y=SD, fill = CellType))+
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "white") +
  scale_fill_manual(values = cellCols)+
  ylab("SD (DNA methylation)")
dev.off()

#----------------------------------------------------------------------#
# POWER CALCULATIONS
#----------------------------------------------------------------------#

plotList <- c( lapply(c(2,5), function(meanDiff){
    allSamples <- calcSamples(allSDs, meanDiff = meanDiff, dataType = "SDs")
    allProps <-calcProps(allSamples)
    x <- plotPower(allProps, "samples")
    return(x)

}),

lapply(c(100, 200), function(sampleSize){
  allSamples <- calcDiff(allSDs, nSamples = sampleSize, dataType = "SDs")
    allProps <-calcProps(allSamples)
    x <- plotPower(allProps, "difference")
    return(x)

}))


combinedPlot <- ggarrange(plotList[[1]], plotList[[2]],plotList[[3]],plotList[[4]], 
                           ncol = 2, nrow = 2, 
                           common.legend = TRUE, legend = "bottom")

pdf(file.path(resPath, paste0("PanelledPowerCurves.pdf"), width = 10, height = 10))
combinedPlot
dev.off()
