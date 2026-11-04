##---------------------------------------------------------------------#
##
## Title: hunk box plot (all celltypes)
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# module load R/4.2.1-foss-2022a


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(tidyr)
library(ggplot2)

dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
outPath <- file.path(dataDir, "4_analysis/plotsForPaper")

cellLabels <- c("NeuN-/SOX10-", "NeuN+", "SOX10+", "Total")
cellTypes <- c("Double-", "NeuN+", "Sox10+", "Total")

setwd(dataDir)


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

# load normalised data
load(normData)

QCmetrics<-QCmetrics[which(QCmetrics$Cell.type %in% cellTypes),] 
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),] 
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename] 


# USCS Annotations DMPs
USCSanno <- read.csv(file.path(outPath, "NeunDMPs9e08UCSCgenes.csv"))
colnames(USCSanno)[1] <- "ProbeID"


#----------------------------------------------------------------------#
# Box plot of top hunk dmp with all cell types
#----------------------------------------------------------------------#

# get probe ID for hunk DMP
hunk <- USCSanno$ProbeID[USCSanno$GeneNames == "HUNK"]

# get P value for hunk DMP
hunkP <- signif(USCSanno$NullModel_SCZ_P[USCSanno$GeneNames == "HUNK"], digits = 2)

plotdf <- as.data.frame(celltypeNormbeta[hunk,])
colnames(plotdf) <- "probe"

if(identical(row.names(plotdf), QCmetrics$Basename)){ # add case control info
    plotdf$Phenotype <- QCmetrics$Phenotype
    plotdf$CellType <- QCmetrics$Cell.type
    } else {
    stop("betas and pheno data don't match!")
    }

plotdf$CellType <- factor(plotdf$CellType, 
                          levels = c("NeuN+", "Sox10+", "Double-", "Total"),
                          labels = c("NeuN+", "SOX10+", "NeuN-/SOX10-", "Total"))

p <- ggplot(plotdf, aes(x = CellType, y=probe, fill=Phenotype)) + 
    geom_boxplot() + 
    geom_point(position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0, dodge.width = 0.75,), alpha = 0.6, size = 1.5) +
    labs(x="Cell Type", y="DNAm (Proportion)", fill = "Phenotype") +
    scale_x_discrete() +
    theme_classic() +
    theme(
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +  
    scale_fill_manual(values = c("Control" = "blue", "Schizophrenia" = "orange")) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
      annotate("text", x = "NeuN+", y = 0.65, label = paste0("P = ", hunkP), fontface = 2) +
         annotate("text", x = "SOX10+", y = 0.8, label = "NS", fontface = 2) +
           annotate("text", x = "NeuN-/SOX10-", y = 0.78, label = "NS", fontface = 2) +
             annotate("text", x = "Total", y = 0.74, label = "NS", fontface = 2)



pdf(file.path(outPath, "boxPlotTopNeuNHunkDMPallCellTypes.pdf"))
print(p)
dev.off()

