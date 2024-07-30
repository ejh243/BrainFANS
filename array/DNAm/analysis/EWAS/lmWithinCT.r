##---------------------------------------------------------------------#
##
## Title: EWAS using linear regression model within each cell type
##
## Purpose of script: perform DNA methylation association analysis of 
## schizophrenia vs controls for each cell type separately and 
## including relevant cell composition covariates
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWASwCC<-function(row,QCmetrics){
   modelNull<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre) 
	modelFull<-lm(row ~ QCmetrics$Phenotype + QCmetrics$Cell.Proportion + QCmetrics$CETYGO + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	modelCC <- lm(row ~ QCmetrics$Phenotype + QCmetrics$Cell.Proportion + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
  # extract case control & cell composition effects
  return(c(
    summary(modelFull)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
    summary(modelFull)$coefficients["QCmetrics$Cell.Proportion",c(1,2,4)],
	summary(modelFull)$coefficients["QCmetrics$CETYGO",c(1,2,4)],
	summary(modelCC)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
    summary(modelCC)$coefficients["QCmetrics$Cell.Proportion",c(1,2,4)],
	summary(modelNull)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)]))	
}

runEWASbasic<-function(row,QCmetrics){
   modelNull<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre) 
   modelFull<-lm(row ~ QCmetrics$Phenotype + QCmetrics$CETYGO + QCmetrics$CETYGO + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	
	  # extract case control & cell composition effects
  return(c(
    summary(modelFull)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
	summary(modelFull)$coefficients["QCmetrics$CETYGO",c(1,2,4)],
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
cellType <- args[2]

normData<-file.path(dataDir, "3_normalised/normalised.rdata")
cellCompData<-file.path(dataDir, "4_analysis", "EstimatedNeuralCellComposition.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)
load(cellCompData)

print(paste0("running EWAS on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]


# add cell composition proportions & CETYGO scores

if (cellType == "Double-"){
   QCmetrics$Cell.Proportion <- predPropBest[QCmetrics$Basename, "NeuNNeg_Sox10Neg_IRF8Pos"]
} else if (cellType == "NeuN+"){
  QCmetrics$Cell.Proportion <- predPropBest[QCmetrics$Basename, "NeuNPos_SOX6Pos"]
}

QCmetrics$CETYGO <- predPropBest[QCmetrics$Basename, "CETYGO"]
 

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]


#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)


if(cellType %in% c("NeuN+", "Double-")){
    clusterExport(cl, list("runEWASwCC"))
	outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWASwCC, QCmetrics), ncol = 18, byrow = TRUE)
	colnames(outtab)<-c(paste0("FullModel", c(rep("_SCZ_", 3), rep("_CellProportion_", 3), rep("_CETYGO_", 3)), c("coeff", "SE", "P")),
		paste0("CCModel", c(rep("_SCZ_", 3), rep("_CellProportion_", 3)), c("coeff", "SE", "P")),
		paste0("NullModel_SCZ_", c("coeff", "SE", "P")))
} else {
    clusterExport(cl, list("runEWASbasic"))
	outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWASbasic, QCmetrics), ncol = 9, byrow = TRUE)
	colnames(outtab)<-c(paste0("FullModel", c(rep("_SCZ_", 3), rep("_CETYGO_", 3)), c("coeff", "SE", "P")),
	paste0("NullModel_SCZ_", c("coeff", "SE", "P")))
}

rownames(outtab)<-rownames(celltypeNormbeta)

save(outtab, file = file.path(resPath, paste0(cellType,"LM.rdata")))
