##---------------------------------------------------------------------#
##
## Title: EWAS using linear regression model within each brain region and cell type
##
## Purpose of script: perform DNA methylation association analysis of 
## schizophrenia vs controls for each cell type and brain region separately
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWASbasicAge<-function(row,QCmetrics){
   modelNull<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$Age + QCmetrics$Sex) 
  return(c(
	summary(modelNull)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)]))	
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1] 
cellType <- args[2] #  IRF8 NeuN SOX10 TN 
brainRegion <- args[3] # Hippocampus Striatum 


normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS")


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)


print(paste0("running EWAS on ", cellType, " cell type and ", brainRegion, " brain region"))

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]

## subset to brain region samples
QCmetrics<-QCmetrics[which(QCmetrics$BrainRegion == brainRegion),]

## rename CaseStatus to Phenotype
QCmetrics$Phenotype <- QCmetrics$CaseStatus

# convert Age to numeric
QCmetrics$Age <- as.numeric(QCmetrics$Age)

# Add in missing age for sample
QCmetrics$Age[is.na(QCmetrics$Age)] <- 71

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]


#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)

print("run EWAS")

clusterExport(cl, list("runEWASbasicAge"))
outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWASbasicAge, QCmetrics), ncol = 3, byrow = TRUE)
colnames(outtab)<-c(paste0("NullModel_SCZ_", c("coeff", "SE", "P")))


rownames(outtab)<-rownames(celltypeNormbeta)

print("saving file")
save(outtab, file = file.path(resPath, paste0(cellType, brainRegion, "LM.rdata")))

