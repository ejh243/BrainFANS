##---------------------------------------------------------------------#
##
## Title: EWAS with mixed effects regression model
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

	pheno<-cbind(row,QCmetrics)
	modelMLM<-lmer(row ~ Phenotype * Cell.type + CCDNAmAge + Sex + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	nullMLM<-lmer(row ~ Phenotype + Cell.type + CCDNAmAge + Sex  +  (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	nullCT<-lmer(row ~ Phenotype + CCDNAmAge + Sex  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	
	# extract case control main effect
	return(c(summary(modelMLM)$coefficients["PhenotypeSchizophrenia",c(1,2,5)],
	
	# extract cell specific case control effect
	anova(modelMLM,nullMLM)[2,8], 
	summary(modelMLM)$coefficients["PhenotypeSchizophrenia:Cell.typeNeuN+",c(1,2,5)],
	summary(modelMLM)$coefficients["PhenotypeSchizophrenia:Cell.typeSox10+",c(1,2,5)],
	
	# extract cell type effect
	anova(nullMLM, nullCT)[2,8],
	summary(nullMLM)$coefficients["Cell.typeNeuN+",c(1,2,5)],
	summary(nullMLM)$coefficients["Cell.typeSox10+",c(1,2,5)]))
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(lme4)
library(lmerTest)
#library(GenABEL)
library(doParallel)
library(devtools)
devtools::load_all(path = "../functionsR")

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
clusterExport(cl, list("runEWAS", "lmer"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 17, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_SCZ_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P", "CellType_P","NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P")  


save(outtab, file = file.path(resPath, "MLM.rdata"))
