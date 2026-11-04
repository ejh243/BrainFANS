##---------------------------------------------------------------------#
##
## Title: Figures for Paper
##
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
library(ggpubr)
library(stringr)
library(ggrepel)
library(qqman)
library(stringi)
library(dplyr)


dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
outPath <- file.path(dataDir, "4_analysis/plotsForPaper")

cellLabels <- c("NeuN-/SOX10-", "NeuN+", "SOX10+", "Total")
cellTypes <- c("Double-", "NeuN+", "Sox10+", "Total")

# Define significance thresholds
pGenomeWide <- 9e-8  # Genome-wide significance
pSuggestive <- 1e-4  # Suggestive significance
FCcutoff <- 0.02  # Fold-change cutoff forMan volcano plots


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

# load results files - annotated with gene info and snps/xhyb probes filtered out
load(file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata"))

# remove 30 probes which have '*' chr and start values of 0/-1
for(i in 1:length(res)){
  res[[i]] <- res[[i]][res[[i]]$chrm != "*", ]   
}

# calculate xlim vals for volcano plot
max_abs_diff <- ceiling(max(sapply(res, function(df) max(abs(df$NullModel_SCZ_coeff), na.rm = TRUE))) * 100) / 100

# read in annototaed DMPs
USCSanno <- read.csv(file.path(outPath, "NeunDMPs9e08UCSCgenes.csv"))
colnames(USCSanno)[1] <- "ProbeID"
USCSanno[USCSanno$GeneNames == ".", "GeneNames"] <- "" # change to blank for no gene



#----------------------------------------------------------------------#
# Specificity scatter plots
#----------------------------------------------------------------------#

# function definition

plotByStatus <- function(dmpCellType, contrastCellType, dmpRes, plotLim, xlabel, ylabel){
    x_col <- paste0("`Coeff:", dmpCellType, "`")
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
dmpRes$sigStatus[which(rownames(dmpRes) %in% USCSanno$ProbeID)] <- "pGenomeWide"


################################## create plots ######################################
######################################################################################

# ensure all plots have same axis sizes
plotLim<-range(dmpRes[,c(4,7,10,13)])

# create plots
p1<- plotByStatus("NeuN+", "Double-", dmpRes, plotLim, "NeuN+", "NeuN-/SOX10-")
p2<- plotByStatus("NeuN+", "Sox10+", dmpRes, plotLim, "NeuN+", "SOX10+")
p3<- plotByStatus("NeuN+", "Total", dmpRes, plotLim, "NeuN+", "Total")

# arrange and save with 9e-08 in red
pdf(file.path(outPath, "ScatterplotsNeuNCelltypeEffectsDiscoveryDMPsLMWithinCTsRedPoints.pdf"), width = 30, height = 10)
ggarrange(p3,p2,p1, nrow = 1, align = "h")
dev.off()




# calculate R values

corDF <-  subset(dmpRes, DiscoveryCellType == "NeuN+")


# calculate R values
cor.test(corDF[, "Coeff:NeuN+"], corDF [, "Coeff:Double-"])
# [1] 0.06218052 , p-value = 0.1992

cor.test(corDF [, "Coeff:NeuN+"], corDF [, "Coeff:Sox10+"])
# [1] 0.2707036  p-value = 1.266e-08

cor.test(corDF [, "Coeff:NeuN+"], corDF [, "Coeff:Total"])
# [1] 0.5441501   p-value < 2.2e-16



################################## violin plots ######################################
######################################################################################

diffLong<-pivot_longer(data.frame(dmpRes[,c(2,4,7,10,13)]), cols = gsub("\\+|-", "\\.", paste("Coeff", cellTypes, sep = ".")))
diffLong <- diffLong %>%
  mutate(
    DiscoveryCellType = stri_replace_all_fixed(DiscoveryCellType, cellTypes, cellLabels, vectorize_all = FALSE),
    name = str_replace_all(name, setNames(cellLabels, unique(diffLong$name))) %>% 
      str_remove("^Coeff\\."))
diffLong$value<-abs(diffLong$value)
colnames(diffLong)[2] <- "CellType"
diffLong$CellType <- factor(diffLong$CellType, levels = c("NeuN+", "SOX10+", "NeuN-/SOX10-", "Total"))
diffLong <- diffLong %>%
  mutate(DiscoveryCellType = factor(DiscoveryCellType, levels = c("NeuN+", "SOX10+", "NeuN-/SOX10-", "Total")))                                        



for (ct in c("NeuN+", "Total")) {
  
  subdat <- diffLong[diffLong$DiscoveryCellType == ct, ]
  
  tiff(file.path(outPath, paste0("ViolinPlot", ct, "DiscoveryDMPs.tiff")), width = 4000, height = 3000, res = 600, compression = "lzw", type='cairo')


  p <- ggplot(subdat, aes(x = CellType, y = value, fill = CellType)) +
    geom_violin() +
    stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
    xlab("") + ylab("Mean difference") +
    #ggtitle(ct) +
    theme_bw() +
    scale_fill_manual(values = cellColors) +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  print(p)
  dev.off()
}


# calculate p values for difference between cell types in each plot
neunStatdat <- dmpRes[dmpRes$DiscoveryCellType == "NeuN+", c("ProbeID", "Coeff:NeuN+", "Coeff:Sox10+", "Coeff:Double-", "Coeff:Total")]

test2 <- wilcox.test(abs(neunStatdat$`Coeff:NeuN+`), abs(neunStatdat$`Coeff:Sox10+`), paired = TRUE, alternative = "greater")
test3 <- wilcox.test(abs(neunStatdat$`Coeff:NeuN+`), abs(neunStatdat$`Coeff:Double-`), paired = TRUE, alternative = "greater")
test4 <- wilcox.test(abs(neunStatdat$`Coeff:NeuN+`), abs(neunStatdat$`Coeff:Total`), paired = TRUE, alternative = "greater")

results <- data.frame(
  comparison        = c("NeuN+ vs SOX10+", "NeuN+ vs Double-", "NeuN+ vs Total"),
  mean_neuron       = round(mean(abs(neunStatdat$`Coeff:NeuN+`)), 4),
  mean_other        = round(c(mean(abs(neunStatdat$`Coeff:Sox10+`)),
                              mean(abs(neunStatdat$`Coeff:Double-`)),
                              mean(abs(neunStatdat$`Coeff:Total`))), 4),
  W_statistic       = c(test2$statistic, test3$statistic, test4$statistic),
  p_value           = c(test2$p.value, test3$p.value, test4$p.value)
)


#----------------------------------------------------------------------#
# Box plots of 16 NeuN+ DMPs
#----------------------------------------------------------------------#


# create plots

plotList <- list()
num=1
for(i in USCSanno$ProbeID){
    plotdf <- as.data.frame(celltypeNormbeta[i,])
    colnames(plotdf) <- "probe"

    if(identical(row.names(plotdf), QCmetrics$Basename)){ # add case control info
        plotdf$Phenotype <- QCmetrics$Phenotype
        } else {
        stop("betas and pheno data don't match!")}

   # create plots
    plotTitle <- paste0(USCSanno[num, "GeneNames"], ": P = ", signif(USCSanno[num,"NullModel_SCZ_P"], 4))
  
    p <- ggplot(plotdf, aes(x = Phenotype, y=probe, fill=Phenotype)) + 
      geom_boxplot() + 
      geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 1.5) +
      labs(title=plotTitle, x=i, y="DNAm (Proportion)")+
      scale_x_discrete() +  # Ensure correct order
      theme_classic() +  # Remove grey background
      theme(
      axis.line = element_line(color = "black"),  # Make axes black
      axis.text = element_text(color = "black"),  # Make axis labels black
      axis.title = element_text(color = "black"), # Make axis titles black
      legend.position = "none") +  
      scale_fill_manual(values = c("Control" = "blue", "Schizophrenia" = "orange")) +
      scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+
      xlab(NULL)

      plotList[[i]] <- p
      num = num + 1
} 


# arrange and save
pdf(file.path(outPath, "boxPlots16NeuNSignificantDMPs.pdf"), width = 20, height = 20)
ggarrange(plotlist = plotList, ncol = 4, nrow = 4)
dev.off()


#----------------------------------------------------------------------#
# Box plot top HUNK DMP on its own
#----------------------------------------------------------------------#

png(file.path(outPath, "boxPlotTopNeuNHunkDMP.png"), width = 2000, height = 2000, res = 300,type = "cairo")
hunk <- plotList[[1]]
hunk <- hunk + theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 20),
                     plot.title = element_text(size = 25, face = "bold")) +
                     ggtitle("")
              
print(hunk)
dev.off()


#----------------------------------------------------------------------#
# VOLCANO PLOTS ALL CELLTYPES
#----------------------------------------------------------------------#

geneNames <- FALSE

# for all cell types
for(i in 1:length(res)){

  volcDat<-as.data.frame(res[[i]])

  volcDat<-volcDat[,c("NullModel_SCZ_coeff","NullModel_SCZ_P","GeneNames")]
  colnames(volcDat)<-c("DNAm_Difference", "P_Value", "Gene")

  # update the significant probes to have the USCS gene names
  matched <- match(rownames(volcDat), USCSanno$ProbeID)
  volcDat$Gene[!is.na(matched)] <- USCSanno$GeneNames[matched[!is.na(matched)]]

# Categorize data into significance groups
volcDat$Differential <- ifelse(volcDat$P_Value < pGenomeWide & abs(volcDat$DNAm_Difference) > FCcutoff, paste0("P<", pGenomeWide),
                          ifelse(volcDat$P_Value < pGenomeWide, paste0("P<", pGenomeWide),
                                 ifelse(volcDat$P_Value < pSuggestive, paste0("P<", pSuggestive),
                                        ifelse(abs(volcDat$DNAm_Difference) > FCcutoff, "BetaDiff", "NS"))))

# Define color scheme
color_scheme <- setNames(
  c("gray", "black", "blue", "red3", "purple"),
  c("NS", "BetaDiff>0.2", "P", paste0("P<", pGenomeWide), paste0("P<", pSuggestive))
)

# Select top genes for annotation
top_genes <- volcDat[volcDat$P_Value < pGenomeWide, ]  # Adjust as needed


# Save plot to tiff
tiff(file.path(outPath, paste0(names(res[i]), "VolcanoPlot.tiff")), width = 4000, height = 3000, res = 600, compression = "lzw", type='cairo')

# Create ggplot volcano plot
p <- ggplot(volcDat, aes(x = DNAm_Difference, y = -log10(P_Value), color = Differential)) +
  geom_point(size = 1.0) +
  scale_color_manual(values = color_scheme) +
  theme_minimal() +
  labs(
    x = "Group Difference in DNAm (Proportion)",
    y = expression(-log[10]("P")),
    color = ""
  ) +
  theme(
    legend.position = "none",
    #legend.position = "right",
    #legend.text = element_text(size = 9),
    panel.grid = element_blank()
  ) +
  # Genome-wide significance threshold
  geom_hline(yintercept = -log10(pGenomeWide), color = "red", linetype = "dashed") +
  # Suggestive significance threshold
  geom_hline(yintercept = -log10(pSuggestive), color = "purple", linetype = "dotted") +
  # Fold change cutoffs
  geom_vline(xintercept = c(-FCcutoff, FCcutoff), color = "black", linetype = "dashed") +
  xlim(max_abs_diff*-1, max_abs_diff) +
  ylim(0,10)

# Add gene annotations if geneNames is TRUE
if(geneNames == TRUE){
xp <- p + geom_label_repel(
  data = top_genes, aes(label = Gene),
  nudge_x = 0.01, nudge_y = 0.4,
  color = "black", size = 2, box.padding = unit(0.2, "lines"), max.overlaps = 70 
)
} else {
  xp <- p
}

# Print plot
print(xp)

dev.off()

}


#----------------------------------------------------------------------#
# MANHATTAN PLOTS FOR ALL CELLTYPES
#----------------------------------------------------------------------#


for(i in 1:length(res)){

  manDat<-as.data.frame(res[[i]])
  manDat <- manDat[,c("chrm","start","NullModel_SCZ_coeff","NullModel_SCZ_P", "GeneNames")]

  colnames(manDat)<-c("CHR", "BP", "DNAm_Difference", "P_Value", "Gene")

  # update the significant probes to have the USCS gene names
  matched <- match(rownames(manDat), USCSanno$ProbeID)
  manDat$Gene[!is.na(matched)] <- USCSanno$GeneNames[matched[!is.na(matched)]]

# Convert CHR to numeric
  manDat$CHR[manDat$CHR == "X"] <- "23" 
  manDat$CHR <- as.numeric(as.character(manDat$CHR))

# Remove any potential NA values in CHR
  manDat <- manDat[!is.na(manDat$CHR), ]
  manDat <- manDat[,c("CHR","BP","P_Value","Gene")]

# Define colnames
  colnames(manDat)<-c("CHR", "BP", "P", "SNP")
  rownames(manDat) <- NULL


# Create Manhattan Plots

  tiff(file.path(outPath, paste0(names(res[i]), "ManhattanPlot.tiff")), width = 4000, height = 3000, res = 600, compression = "lzw", type='cairo')

  par(mfrow = c(1,1))
  par(mar = c(4, 4, 2, 1))
  manhattan(manDat,
          col = c("darkblue", "lightblue"), # Alternating colors
          genomewideline = -log10(pGenomeWide), # Genome-wide significance threshold
          suggestiveline = -log10(pSuggestive), # Suggestive threshold
          annotatePval = pGenomeWide,
          annotateTop = FALSE,
          cex.axis = 0.6,
          ylim = c(0, 10))
dev.off()

}


#----------------------------------------------------------------------#
# QQ PLOTS ALL CELLTYPES
#----------------------------------------------------------------------#

## ggplot with sig DMPs in red

cellTypeOrder  <- c("NeuN+", "Sox10+", "Double-", "Total")
cellTypeLabels <- c("NeuN+", "SOX10+", "NeuN-/SOX10-", "Total")


plotList <- list()

for(i in 1:length(cellTypeOrder)){

  # calculate lambda
    qqDat <- as.data.frame(res[[cellTypeOrder[i]]])
    qqDat$sig <- ifelse(qqDat$NullModel_SCZ_P < pGenomeWide, "Sig", "NS")
    qqDat <- qqDat[order(qqDat$NullModel_SCZ_P), ]

    qqDat$expected <- -log10(ppoints(length(qqDat$NullModel_SCZ_P)))
    qqDat$observed <- -log10(qqDat$NullModel_SCZ_P)

    alpha <- median(qchisq(1-qqDat$NullModel_SCZ_P, 1))/qchisq(0.5,1)


# plot
plotList[[i]] <- ggplot(qqDat, aes(x = expected, y = observed)) +
  geom_point(aes(colour = sig)) +
  scale_colour_manual(
    values = c("NS" = "black", "Sig" = "red")
  ) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 12), expand = FALSE) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)#,
   # plot.margin = unit(c(1,1,1,2), "cm")
    ) +
  annotate(
    "text",
    x = 2,
    y = 10,
    label = paste("lambda =", signif(alpha, 3)),
    size = 4
  ) +
  ggtitle(cellTypeLabels[i])

}



# save as one plot with shared axis labels
figure <- ggarrange(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[4]], 
                    ncol = 2, nrow = 2, align = "h")

figure <- annotate_figure(figure,
  left = text_grob("Observed -log10(p-value)", rot = 90, size = 14),
  bottom = text_grob("Expected -log10(p-value)", size = 14))


png(file.path(outPath, "QQPlotsAllCellTypes_ggplot.png"), width = 2000, height = 2000, res = 300,type = "cairo")
figure
dev.off()


