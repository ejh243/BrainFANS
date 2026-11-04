##---------------------------------------------------------------------#
##
## Title: All Brain Region BoxPlots
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# ignores epicV2 duplicate probes in the Hippocampus and Striatum data
# only keeps probes that are in both epicv1 and epicv2 datasets


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)

PFCcellType <- "NeuN+"
HScellType <- "NeuN"

dataDir <- "Schizophrenia/brainRegionData/"
outPath <- file.path(dataDir, "4_analysis/plots")

HSpath <- file.path(dataDir, "3_normalised/normalised.rdata")
newPFCpath <- file.path(dataDir, "3_normalised/PFCnormalised.rdata")

v2manifest <- "references/EPICArray/EPICv2_reannotated_manifest_v1.0.csv"

allRes <- file.path(outPath, paste0(PFCcellType,"FullResultsAllBrainRegions.csv"))


#----------------------------------------------------------------------#
# LOAD AND PREPARE PFC BETAS MATRIX
#----------------------------------------------------------------------#

load(newPFCpath)

# subset to just correct cell type and samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == PFCcellType),]
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]
PFCQCmetrics <- QCmetrics
PFCQCmetrics$BrainRegion <- "PFC"
rm(QCmetrics)

PFCbeta<-celltypeNormbeta[,PFCQCmetrics$Basename]
rm(celltypeNormbeta)

PFCbeta <- as.data.frame(PFCbeta)
PFCbeta$ProbeID <- row.names(PFCbeta)


#----------------------------------------------------------------------#
# LOAD AND PREPARE HIPPOCAMPUS AND STRIATUM BETAS MATRIX
#----------------------------------------------------------------------#

load(HSpath)

# subset to just correct cell type and samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == HScellType),]
QCmetrics<-QCmetrics[QCmetrics$CaseStatus %in% c("Schizophrenia", "Control"),]
HSQCmetrics <- QCmetrics
colnames(HSQCmetrics)[6] <- "Phenotype"
rm(QCmetrics)

HSbeta<-celltypeNormbeta[,HSQCmetrics$Basename]
rm(celltypeNormbeta)

HSbeta <- as.data.frame(HSbeta)
HSbeta$IlmnID <- row.names(HSbeta)

# load epicv2 manifest to match IlmnID to EPICv1_Loci
manifest <- read.csv(v2manifest, stringsAsFactors = FALSE)
manifest <- manifest[, c("IlmnID", "EPICv1_Loci")]
colnames(manifest)[2] <- "ProbeID"

# join with manifest to get EPICv1_Loci
HSbeta <- dplyr::left_join(HSbeta, manifest, by = "IlmnID")
HSbeta<- HSbeta[!duplicated(HSbeta$ProbeID),] # remove duplicates
row.names(HSbeta) <- HSbeta$ProbeID # set row names to ProbeID for consistency


#----------------------------------------------------------------------#
# Boxplots of NeuN DMPs
#----------------------------------------------------------------------#

# load results object
res <- read.csv(allRes, stringsAsFactors = FALSE)
row.names(res) <- res$ProbeID # set row names to ProbeID for consistency

# extact unique genes
res$unique_genes <- sapply(res$GeneNames, function(x) {
unique_genes <- unique(unlist(strsplit(x, ";")))
        paste(unique_genes, collapse = ";")            
            })  

# manually add gene name for SSC5D probe
res["cg12421899", "unique_genes"] <- "SSC5D"            



# Plot

brainRegions <- c("fullPFC", "Hipp", "Striatum")

for(i in brainRegions){
# Get top 20 DMPs for region
    DMPs <- res[order(res[[paste0("P.", i)]]), ][1:20, "ProbeID"]

    PFCplotBetas <- PFCbeta[DMPs, -ncol(PFCbeta)] # subset betas matrix to top 20 probes and remove probeID col
    HSplotBetas <- HSbeta[DMPs, 1:(ncol(HSbeta)-2)] # subset betas matrix to top 20 probes

    plotBetas <- cbind(PFCplotBetas, HSplotBetas)
    plotQCmetrics <- rbind(PFCQCmetrics[, c("Basename", "Phenotype", "BrainRegion")], HSQCmetrics[, c("Basename", "Phenotype", "BrainRegion")])

    identical(colnames(plotBetas), plotQCmetrics$Basename) # check betas and QC metrics match


    ## make plots
    
    plotList <- list()
    for(j in DMPs){
    
        plotdf <- as.data.frame(t(plotBetas[j,]))
        colnames(plotdf) <- "probe"
        #plotdf$probe <- as.numeric(plotdf$probe) # convert to numeric

        if(identical(row.names(plotdf), plotQCmetrics$Basename)){ # add pheno info
         plotdf <- cbind(plotdf, plotQCmetrics)
            } else {
         stop("betas and pheno data don't match!")
            }
        plotdf$BrainRegion <- factor(plotdf$BrainRegion, levels = c("PFC", "Hippocampus", "Striatum"))

    # create plot title
        plotTitle <- res[j, "unique_genes"]

    
        p <- ggplot(plotdf, aes(x = BrainRegion, y=probe, fill=Phenotype)) + 
        geom_boxplot() + 
    geom_point(position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0, dodge.width = 0.75,), alpha = 0.6, size = 1.5) +
        labs(title=plotTitle, x=NULL, y="DNAm (Proportion)")+
        scale_x_discrete() +  # Ensure correct order
  theme_classic() +  # Remove grey background
  theme(
    axis.line = element_line(color = "black"),  # Make axes black
    axis.text = element_text(color = "black"),  # Make axis labels black
    axis.title = element_text(color = "black")) +  # Make axis titles black
     scale_fill_manual(values = c("Control" = "blue", "Schizophrenia" = "orange"))+
     theme(plot.title = element_text(hjust = 0.5))
       #theme(legend.position = "none") # Remove legend
  

plotList[[j]] <- p

    }

    png(file.path(outPath, paste0(PFCcellType, "matched", i, "top2DMPsBoxplots.png")), width = 3000, height = 1500, res = 300,type = "cairo")
    print(ggarrange(plotList[[1]], plotList[[2]], ncol = 2, nrow = 1, common.legend = TRUE, legend="right"))    
    dev.off()
}


