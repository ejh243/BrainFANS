##---------------------------------------------------------------------#
##
## Title: GWAS/SCHEMA GENES DMPs Boxplots
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS AND LOAD DATA
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)

# function to extract unique genes
uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }


dataDir <- "Schizophrenia"
refDir <- "references/"
normData <- paste0(dataDir, "/3_normalised/normalised.rdata")


setwd(dataDir)

load(normData)

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to cell types used
cellTypes <- c("NeuN+")
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]

# subset betas matrix to same samples
celltypeNormbeta <- celltypeNormbeta[, QCmetrics$Basename]


# load EWAS results 
load(file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata"))
outtab <- res[[2]]


# load gene enrichment results
dmps <- read.csv(file.path(dataDir, "geneEnrichment_NeuN+LM_1e04HG38_DMPs_Combined.csv"))
cpgs <- dmps$IlmnID



# extract unique genes
outtab$unique_genes <- sapply(outtab$GeneNames, function(x) {
  unique_genes <- unique(unlist(strsplit(x, ";")))
  paste(unique_genes, collapse = ";")            
})  



#----------------------------------------------------------------------#
# Boxplots of NeuN DMPs
#----------------------------------------------------------------------#

plotBetas <- celltypeNormbeta[cpgs, QCmetrics$Basename] # subset betas matrix 


plotList <- list()
for(i in cpgs){

    plotdf <- as.data.frame(plotBetas[i,])
    colnames(plotdf) <- "probe"

    if(identical(row.names(plotdf), QCmetrics$Basename)){ # add pheno info
                    plotdf$Phenotype <- QCmetrics$Phenotype
                    } else {
                        stop("betas and pheno data don't match!")
            }

    
    p <- ggplot(plotdf, aes(x = Phenotype, y=probe, fill=Phenotype)) + 
      geom_boxplot() + 
      geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.5) +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 10)) +
      labs(x="", y="DNAm (Proportion)")+
      scale_x_discrete() +  # Ensure correct order
      theme_classic() +  # Remove grey background
      theme(legend.position = "none",
            axis.line = element_line(color = "black"),  # Make axes black
            axis.text = element_text(color = "black", size = 12),  # Make axis labels black
            axis.title = element_text(color = "black", size = 14)
            ) +  # Make axis titles black
      scale_fill_manual(values = c("Control" = "blue", "Schizophrenia" = "orange"))


    plotList[[i]] <- p

} 


png(file.path(dataDir, "4_analysis/geneSetEnrichment/NeunDMPGWASSCHEMABoxplots.png"), width = 3000, height = 3000, res = 300,type = "cairo")
ggarrange(plotlist = plotList, ncol = 3, nrow = 3)
dev.off()


#----------------------------------------------------------------------#
# Output DMPs
#----------------------------------------------------------------------#

colsToKeep <- c("IlmnID", "NullModel_SCZ_coeff", "NullModel_SCZ_SE", "NullModel_SCZ_P", "CHR", "MAPINFO",                  
"Gene", "sczGWAS.Gene", "SCHEMA.Gene")

dmps <- dmps[,colsToKeep]

colnames(dmps)[2:4] <- c("Effect_Size", "SE", "P")

write.csv(dmps, file.path(dataDir, "4_analysis/geneSetEnrichment/GWASandSchemaDMPList.csv"))


