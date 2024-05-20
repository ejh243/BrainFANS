##---------------------------------------------------------------------#
##
## Title: EWAS with linear regression model
##
## Purpose of script: perform DNA methylation association analysis of 
## schizophrenia vs controls testing for main and cell-specific effects.
## simulataneously testing for cell type differences
##
## Author: Eilis Hannon
##
## Date Created: 2022-07-08
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS<-function(row,QCmetrics){

	modelLM<-lm(row ~ QCmetrics$Phenotype * QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	nullLM<-lm(row ~ QCmetrics$Phenotype + QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	nullCT<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	
	# extract case control main effect
	return(c(summary(modelLM)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
	
	# extract cell specific case control effect
	anova(modelLM,nullLM)[2,6], 
	summary(modelLM)$coefficients["QCmetrics$PhenotypeSchizophrenia:QCmetrics$Cell.typeNeuN+",c(1,2,4)],
	summary(modelLM)$coefficients["QCmetrics$PhenotypeSchizophrenia:QCmetrics$Cell.typeSox10+",c(1,2,4)],
	
	# extract cell type effect
	anova(nullLM, nullCT)[2,6],
	summary(nullLM)$coefficients["QCmetrics$Cell.typeNeuN+",c(1,2,4)],
	summary(nullLM)$coefficients["QCmetrics$Cell.typeSox10+",c(1,2,4)]))
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

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis/EWAS")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
cellTypes<-unique(QCmetrics$Cell.type)

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 17, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_SCZ_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P", "CellType_P","NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P")  


save(outtab, file = file.path(resPath, "LM.rdata"))