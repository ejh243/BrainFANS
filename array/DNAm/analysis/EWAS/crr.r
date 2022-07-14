##---------------------------------------------------------------------#
##
## Title: EWAS with clustered robust regression model
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

	p.df <- pdata.frame(data.frame("meth" = row, "phenotype" = QCmetrics$Phenotype, "age" = QCmetrics$CCDNAmAge, "sex" = QCmetrics$Sex, "cell.type" = QCmetrics$Cell.type, "brain.bank" = QCmetrics$Tissue.Centre, "id" = QCmetrics$Indidivual.ID), index = c("id"), drop.index = F)
	
	modelP <- plm(meth ~ phenotype * cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	nullP <- plm(meth ~ phenotype + cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	nullCT<-plm(meth ~ phenotype + age + sex + brain.bank, data = p.df, model = "pooling")
	
	# compute Stata like df-adjustment
	G <- length(unique(p.df$id))
	N <- length(p.df$id)
	dfa <- (G/(G - 1)) * (N - 1)/modelP$df.residual
	
	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(modelP, type = "HC0", cluster = "group", adjust = T)
	model<-coeftest(modelP, vcov = firm_c_vcov)

	# extract case control main effect
	c(model["phenotypeSchizophrenia",c(1,2,4)],
	
	# extract cell specific case control effect
	waldtest(modelP,nullP, vcov = firm_c_vcov, test = "F")[2,4],
	model["phenotypeSchizophrenia:cell.typeNeuN+",c(1,2,4)],
	model["phenotypeSchizophrenia:cell.typeSox10+",c(1,2,4)],
	
	# extract cell type effect
	waldtest(modelP, nullCT, vcov = firm_c_vcov, test = "F")[2,4],
	model["cell.typeNeuN+",c(1,2,4)],
	model["cell.typeSox10+",c(1,2,4)])	

}




#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(plm)
library(lmtest)
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
clusterExport(cl, list("runEWAS", "lmer", "pdata.frame", "plm", "vcovHC", "coeftest", "waldtest"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 17, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_SCZ_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P", "CellType_P","NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P")  


save(outtab, file = file.path(resPath, "RobustClusteredErrors.rdata"))