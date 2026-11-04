
##---------------------------------------------------------------------#
##
## Title: PFC vs Other Brain Region Plots and tables
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
library(dplyr)
library(ggpubr)

dataDir <- "Schizophrenia/brainRegionData/"
resDir <- file.path(dataDir, "4_analysis/EWAS")
outPath <- file.path(dataDir, "4_analysis/plots")

v2manifest <- "references/EPICArray/EPICv2_reannotated_manifest_v1.0.csv"

hippPath <- file.path(resDir, "NeuNHippocampusLM.rdata")
striatumPath <- file.path(resDir, "NeuNStriatumLM.rdata")
newPFCPath <- file.path(resDir, "NeuN+PFCLM.rdata")
PFCpath <- file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata")

# define thresholds
pGenomeWide <- 9e-8  # Genome-wide significance
pSuggestive <- 1e-4  # Suggestive significance


#----------------------------------------------------------------------#
# LOAD AND PREPARE Striatum and Hippocampus Data
#----------------------------------------------------------------------#

# load epicv2 manifest to match IlmnID to EPICv1_Loci
manifest <- read.csv(v2manifest, stringsAsFactors = FALSE)
manifest <- manifest[, c("IlmnID", "EPICv1_Loci")]
colnames(manifest)[2] <- "ProbeID"

# load hippocampus, striatum and new PFC  results and combine
load(hippPath)
hipp <- as.data.frame(outtab)
rm(outtab)

load(striatumPath)
striatum <- as.data.frame(outtab)
rm(outtab)

regionRes <- cbind(hipp, striatum)
regionRes$IlmnID <- rownames(regionRes)

# reorder columns for combining with PFC res
matchData <- regionRes[, c(3,1,2,6,4,5,7)]
colnames(matchData) <- c("P:Hipp", "Coeff:Hipp", "SE:Hipp",
                         "P:Striatum", "Coeff:Striatum", "SE:Striatum",
                         "IlmnID")

# join with manifest to get EPICv1_Loci
matchData <- dplyr::left_join(matchData, manifest, by = "IlmnID")
matchData <- matchData[!duplicated(matchData$ProbeID),] # remove duplicates

# load new PFC EWAS results
load(newPFCPath)
newPFC <- as.data.frame(outtab)
rm(outtab)

newPFC$ProbeID <- rownames(newPFC)
newPFC <- newPFC[, c("ProbeID", "NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")]
colnames(newPFC) <- c("ProbeID", "P:newPFC", "Coeff:newPFC", "SE:newPFC")


#----------------------------------------------------------------------#
# LOAD AND PREPARE PFC Data
#----------------------------------------------------------------------#

# load PFC results files - annotated with gene info and snps/xhyb probes filtered out
load(PFCpath)

# remove 30 probes which have '*' chr and start values of 0/-1
for(i in 1:length(res)){
  res[[i]] <- res[[i]][res[[i]]$chrm != "*", ]   
}

cellTypes <- c("Double-", "NeuN+", "Sox10+", "Total")

# load NEUN DMPs - list of probes significant at pGenomeWide threshold
USCSanno <- read.csv(file.path(dataDir, "4_analysis/plotsForPaper/NeunDMPs9e08UCSCgenes.csv"))
colnames(USCSanno)[1] <- "ProbeID"
USCSanno[USCSanno$GeneNames == ".", "GeneNames"] <- USCSanno[USCSanno$GeneNames == ".", "ProbeID"]


# create single dataframe of all probes significant at threshold for all the cell types
thres<-1e-04 

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

dim(dmpRes)



#----------------------------------------------------------------------#
# Join PFC dmp results with new PFC res, hippocampus and striatum results
#----------------------------------------------------------------------#

dmpRes <- dplyr::left_join(dmpRes, matchData, by = "ProbeID")
dmpRes <- dplyr::left_join(dmpRes, newPFC, by = "ProbeID")

dmpRes <- dmpRes[complete.cases(dmpRes), ] # remove rows with NAs

# add extra colum for colouring significant DMPs
dmpRes$sigStatus <- rep("pSuggestive")
dmpRes$sigStatus[which(dmpRes[,"ProbeID"] %in% USCSanno$ProbeID)] <- "pGenomeWide"

# remove IlmnID column
dmpRes <- dmpRes[, !(names(dmpRes) %in% "IlmnID")]


#----------------------------------------------------------------------#
# Specificity scatter plots
#----------------------------------------------------------------------#

# function definition

plotByStatus <- function(dmpCellType, Xcol, contrastCellType, dmpRes, plotLim, xlabel, ylabel){
    x_col <- paste0("`Coeff:", Xcol, "`")
    y_col <- paste0("`Coeff:", contrastCellType, "`")

    subset(dmpRes, DiscoveryCellType == dmpCellType) %>%
    ggplot(aes_string(x = x_col, y = y_col, colour = "sigStatus")) + 
    geom_point() +
    xlab(xlabel) + ylab(ylabel) + 
    #ggtitle(xlabel) +
    xlim(plotLim) + ylim(plotLim) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +

theme_minimal() +
  theme(
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 35, face = "bold"),
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


#----------------------------------------------------------------------#
# Create and Save Plots
#----------------------------------------------------------------------#

# for 1e-04 threshold 
dmpRes <- dmpRes[which(dmpRes[,6] < 1e-04),]

# ensure all plots have same axis sizes
plotLim <- range(dmpRes[,c(4,7,10,13,16,19,22)])
  
# create plots
p1 <- plotByStatus("NeuN+", "newPFC", "Hipp", dmpRes, plotLim, "PFC", "Hippocampus") + annotate("text", x=-0.15, y=0.11, label= "R = 0.41", size=14)
p2 <- plotByStatus("NeuN+", "newPFC", "Striatum", dmpRes, plotLim, "PFC", "Striatum") + annotate("text", x=-0.15, y=0.11, label= "R = 0.19", size=14)
p3 <- plotByStatus("NeuN+", "newPFC", "NeuN+", dmpRes, plotLim, "PFC (matched samples)", "PFC (all samples)") + annotate("text", x=-0.15, y=0.11, label= "R = 0.95", size=14)
  
# arrange and save
pdf(file.path(outPath, "1e-04ScatterplotsNeuN_BrainRegionsRedPoints.pdf"), width = 30, height = 10)
print(ggarrange(p3, p1, p2, nrow = 1, align = "h"))

dev.off()



