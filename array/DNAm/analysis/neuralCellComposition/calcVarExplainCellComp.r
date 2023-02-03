##---------------------------------------------------------------------#
##
## Title: Variance explained by neural cell composition 
##
## Purpose of script: To calculate variance explained by cellular 
## composition, and other biological and technical factors.
##
## Author: Eilis Hannon
##
## Date Created: 13/12/2022
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# assumes cellular compoisition already calculated
# data provided as R object with betas matrix and pheno data.frame
# excludes samples < 18


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)

bulkPath <- args[1]
cellPath <- args[2]
resPath <- args[3]
outName<- args[4]

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(variancePartition)


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#
load(bulkPath)
load(cellPath)

# apply to PCAs
pcaBetas <- prcomp(na.omit(t(betas[, which(pheno$Age > 17)])))
pcaBetasProp <- pcaBetas$sdev^2/sum(pcaBetas$sdev^2)
pcaBetasScore <- pcaBetas$x[,which(pcaBetasProp > 0.01)]

pcaPartIDOL<-list()
pcaPartANOVA<-list()
varPartIDOL<-list()
varPartANOVA<-list()
for(i in 1:length(predCCIDOL)){
	cellTypes <- colnames(predCCANOVA[[i]])[1:(ncol(predCCANOVA[[i]])-2)]
	form <- as.formula(paste("~ Age + Sex + Sentrix_ID + CETYGO ", paste(cellTypes, collapse = " + "),sep = " + "))
	pheno.tmp<-cbind(pheno, predCCANOVA[[i]][pheno$Sentrix_Full,c(cellTypes, "CETYGO")])
	pheno.tmp<-pheno.tmp[which(pheno.tmp$Age > 17),]
	pcaPartANOVA[[i]] <- fitExtractVarPartModel(t(pcaBetasScore), form, pheno.tmp)
	varPartANOVA[[i]] <- fitExtractVarPartModel(betas, form, pheno.tmp)

	if(!is.null(predCCIDOL[[i]])){
		pheno.tmp<-cbind(pheno, predCCANOVA[[i]][pheno$Sentrix_Full,c(cellTypes, "CETYGO")])
		pcaPartIDOL[[i]] <- fitExtractVarPartModel(t(pcaBetasScore), form, pheno.tmp)
		varPartIDOL[[i]] <- fitExtractVarPartModel(betas, form, pheno.tmp)
	}
}                              

save(varPartIDOL, varPartANOVA, pcaPartIDOL, pcaPartANOVA, file = file.path(resPath, "varExplained", paste0(outName, "VarExplainedCellComp.Rdata")))

