##---------------------------------------------------------------------#
##
## Title: Null vs Full Model Coefficients and CC Model Coefficients
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)


dataDir <- "Schizophrenia"
outPath <- file.path(dataDir, "4_analysis/plotsForPaper")

# full results, with snps etc removed
load(file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata"))

# cell types
cellTypes <- "NeuN+"

## discovery threshold
thres<-1e-4


#----------------------------------------------------------------------#
# IDENTIFY DMPS
#----------------------------------------------------------------------#

# just neun
neun <- res[[cellTypes]]

# discovery dmps
neun <- neun[which(neun$NullModel_SCZ_P < thres),]


#----------------------------------------------------------------------#
# Plot CC vs null model coefficients for discovery dmps in NeuN+
#----------------------------------------------------------------------#


pCC <- ggplot(neun, aes(x = NullModel_SCZ_coeff, y = CCModel_SCZ_coeff)) +
    geom_point() +
    xlab("NullModel") + 
    ylab("CCModel") + 
    ggtitle("NeuN+") +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +

    theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none") + 
geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40")


png(file.path(outPath, "CCModelvsNullModel_NeuN+.png"), width = 2000, height = 2000, res = 300,type = "cairo")
pCC
dev.off()


# R value
cor(neun$NullModel_SCZ_coeff, neun$CCModel_SCZ_coeff)
# 0.9958224

# P value
neun$Nullsign <- sign(neun$NullModel_SCZ_coeff)
neun$CCsign <- sign(neun$CCModel_SCZ_coeff)

CCNullSignmatch <- neun$Nullsign == neun$CCsign

  
signtest <- binom.test(sum(CCNullSignmatch == T), length(CCNullSignmatch), p = 0.5, alternative = "greater")
signtest$p.value
# 1.442653e-129
