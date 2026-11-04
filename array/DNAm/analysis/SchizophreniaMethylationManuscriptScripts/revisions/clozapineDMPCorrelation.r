## ---------------------------------------------------------------------#
##
## Title: Correlate DMPs with clozapine results
##
## ---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggrepel)
library(qqman)
library(stringi)


dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
outPath <- file.path(dataDir, "4_analysis/clozapineOverlaps")

cellLabels <- c("NeuN-/SOX10-", "NeuN+", "SOX10+", "Total")
cellTypes <- c("Double-", "NeuN+", "Sox10+", "Total")

# Define significance thresholds
pGenomeWide <- 9e-8  # Genome-wide significance
pSuggestive <- 1e-4  # Suggestive significance


# define cell colours
cellColors <- c(
  "NeuN-/SOX10-" = "#E69F00",
  "NeuN+"        = "#56B4E9",
  "SOX10+"       = "#009E73",
  "Total"        = "#CC79A7"
)

setwd(dataDir)


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

# load normalised data
load(normData)

QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == "NeuN+"),] ## subset to NeuN+
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),] # filter to schizophrenia and control only
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename] # subset beta matrix to analysis samples

# load results files - annotated wiht gene info and snps/xhyb probes filtered out
load(file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata"))

# remove 30 probes which have '*' chr and start values of 0/-1
for(i in 1:length(res)){
  res[[i]] <- res[[i]][res[[i]]$chrm != "*", ]   
}


# load clozapine results
fullClozRes <- read.csv("Clozapine/ResultsTablesForPaper/MergedFullEWASResults.csv", row.names = 1)


#----------------------------------------------------------------------#
# Scatter plots
#----------------------------------------------------------------------#

# function definition

plotByStatus <- function(dmpCellType, clozDat, dmpRes, plotLim, xlabel, ylabel){

  # subset to cell type
  dmpSub <- subset(dmpRes, DiscoveryCellType == dmpCellType)
  
  # merge with clozapine results on probe ID
  plotDat <- merge(dmpSub, clozDat, by = "ProbeID") 
  print(paste0("N overlap DMPS: ", nrow(plotDat)))
  
  x_col <- paste0("Coeff:", dmpCellType)
  y_col <- "Linear_EffectSize"
  
  sameDir <- sum(sign(plotDat$`Coeff:NeuN+`) == sign(plotDat$Linear_EffectSize))
  Pval <- binom.test(sameDir, nrow(plotDat), p = 0.5)

  Rval <- cor.test(plotDat[[x_col]], plotDat[[y_col]])

  
  
  ggplot(plotDat, aes_string(x = paste0("`", x_col, "`"), y = y_col, colour = "sigStatus")) + 
    geom_point() +
    xlab(xlabel) + ylab(ylabel) + 
    xlim(plotLim) + ylim(plotLim) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    ggtitle(paste0("r = ", signif(Rval$estimate, 3), "   p = ", signif(Pval$p.value, 3))) +
    
    theme_minimal() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 15, face = "bold"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.6),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c(
      "pGenomeWide" = "red",
      "pSuggestive" = "black"))          
}



# create single dataframe of all probes significant at threshold for all the cell types
thres<-pSuggestive

dmpList<-NULL
for(i in 1:4){
  dmpList<-rbind(dmpList, cbind(rownames(res[[i]])[which(res[[i]][,"NullModel_SCZ_P"] < thres)], cellTypes[i]))
}
colnames(dmpList)<-c("ProbeID", "DiscoveryCellType")


dmpRes<-cbind(dmpList, res[[1]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")],
              res[[2]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")],
              res[[3]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")],
              res[[4]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")])
colnames(dmpRes)<-c("ProbeID", "DiscoveryCellType", outer(c("P", "Coeff", "SE"), cellTypes, paste, sep = ":"))

# add extra colum for colouring significant DMPs
dmpRes$sigStatus <- rep("pSuggestive")
dmpRes$sigStatus[which(dmpRes$`P:NeuN+` < 9e-08)] <- "pGenomeWide"




################################## create plots ######################################
######################################################################################

# ensure all plots have same axis sizes
plotLim <- c(-0.16, 0.16)

# merged res
p1 <- plotByStatus("NeuN+", fullClozRes, dmpRes, plotLim, "NeuN+", "Clozapine (450K/EPIC merge)")


# save
pdf(file.path(outPath, "discoveryDMPmergedClozapineScatterplot.pdf"))
print(p1)
dev.off()


