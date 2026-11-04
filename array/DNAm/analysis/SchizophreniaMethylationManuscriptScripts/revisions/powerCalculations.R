##---------------------------------------------------------------------#
##
## Title: Power calculations
##
##
## Input files: normalised betas and QCmetrics 
##
##
##---------------------------------------------------------------------#



#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

gc()


library(CellPower)

dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")

outDir <- file.path(dataDir, "4_analysis/powerCalcs")


cellTypes <- c("Double-", "NeuN+", "Sox10+", "Total")


#----------------------------------------------------------------------#
# LOAD AND FORMAT DATA
#----------------------------------------------------------------------#

load(normData)

# filter 
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,] ## subset to cell types used
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),] # filter to schizophrenia and control only
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename] # subset beta matrix to analysis samples


#----------------------------------------------------------------------#
# POWER - n samples
#----------------------------------------------------------------------#

# to calculate the number of samples needed to detect a specified difference
allSamples <- calcSamples(celltypeNormbeta, QCmetrics, meanDiff = 5, dataType = "betaMatrix")

# to calculate the proportion of probes with power > 0.8
allProps <-calcProps(allSamples)

# to plot the power curve 
myPlot <- plotPower(allProps, "samples")

# save 
pdf(file.path(outDir, "samplesPlot.pdf"))
print(myPlot)
dev.off()






