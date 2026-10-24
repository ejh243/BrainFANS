##---------------------------------------------------------------------#
##
## Title: EWAS with mixed effects regression model
##
## Purpose of script: perform DNA methylation association analysis of 
## schizophrenia vs controls testing for main and cell-specific effects.
## simulataneously testing for cell type differences
##
##
##---------------------------------------------------------------------#

# module load R/4.2.1-foss-2022a

# module load Anaconda3
# source activate EWAS_parallel

#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS<-function(row,QCmetrics){

	pheno<-cbind(row,QCmetrics)
	modelMLM<-lmer(row ~ Phenotype * Cell.type + Age + Sex + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	nullMLM<-lmer(row ~ Phenotype + Cell.type + Age + Sex  +  (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	nullCT<-lmer(row ~ Phenotype + Age + Sex  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	
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
library(doParallel)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- "Schizophrenia"

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis/EWAS/interactionModel")

pSuggestive <- 1e-4  # Suggestive significance


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)


# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to cell types used
cellTypes <- c("Double-", "NeuN+", "Sox10+")
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]


celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename] # subset beta matrix to analysis samples
betas <- celltypeNormbeta


#----------------------------------------------------------------------#
# GET DMPs
#----------------------------------------------------------------------#

# load results files - annotated with gene info and snps/xhyb probes filtered out
load(file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata"))

# remove 30 probes which have '*' chr and start values of 0/-1
for(i in 1:length(res)){
  res[[i]] <- res[[i]][res[[i]]$chrm != "*", ]   
}

# create single dataframe of all probes significant at threshold for all the cell types
thres<-pSuggestive

dmpList<-NULL
for(i in 1:3){
    dmpList<-rbind(dmpList, cbind(rownames(res[[i]])[which(res[[i]][,"NullModel_SCZ_P"] < thres)], cellTypes[i]))
}
colnames(dmpList)<-c("ProbeID", "DiscoveryCellType")


dmpRes<-cbind(dmpList, res[[1]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")],
res[[2]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")],
res[[3]][dmpList[,1], c("NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")])
colnames(dmpRes)<-c("ProbeID", "DiscoveryCellType", outer(c("P", "Coeff", "SE"), cellTypes, paste, sep = ":"))


#----------------------------------------------------------------------#
# NEUN DMPS ONLY
#----------------------------------------------------------------------#

nOnly <- dmpRes$ProbeID[dmpRes$DiscoveryCellType == "NeuN+"]
celltypeNormbeta <- betas[nOnly,] # 428

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "lmer"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 17, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_SCZ_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P", "CellType_P","NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P")  


save(outtab, file = file.path(resPath, "neunOnlyMLM.rdata"))


# discovery threshold
outtab <- as.data.frame(outtab)
sum(outtab$CellType_SCZ_P < 0.05)
max(outtab$CellType_SCZ_P)
max(outtab$CellType_P)


# genome wide threshold 
n16 <- dmpRes$ProbeID[which(dmpRes$`P:NeuN+` < 9e-08)] 
nRes <- outtab[n16,]

sum(nRes$CellType_SCZ_P < 0.05) 
max(nRes$CellType_SCZ_P) 
max(nRes$CellType_P)


#----------------------------------------------------------------------#
# ALL CELLTYPE DMPS
#----------------------------------------------------------------------#

all <- dmpRes$ProbeID
celltypeNormbeta <- betas[all,] 

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 17, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_SCZ_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P", "CellType_P","NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P")  


save(outtab, file = file.path(resPath, "allCellTypes1e-4_MLM.rdata"))


outtab <- as.data.frame(outtab)
sum(outtab$CellType_SCZ_P < 0.05)
max(outtab$CellType_SCZ_P) 
max(outtab$CellType_P)

