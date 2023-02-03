##---------------------------------------------------------------------#
##
## Title: Plot variance explained statistics 
##
## Purpose of script: To summarise variance explained by cellular 
## composition.
##
## Author: Eilis Hannon
##
## Date Created: 03/02/2023
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
## set plotting colours
library(paletteer)
ctCols <- paletteer_d("ggsci::category10_d3")
names(ctCols)<-c("DoubleNeg", "NeuNPos", "NEUNNeg", "Sox10Pos", "IRF8Pos", "TripleNeg", "SATB2Neg","SATB2Pos", 
"SOX6Neg", "SOX6Pos")


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)

resPath <- args[1]
outName<- args[2]

load(file.path(resPath, paste0(outName, "VarExplainedCellComp.Rdata"))

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(variancePartition)

# caluclate total variance explained by all cell types
for(i in 1:length(varPartANOVA)){
    varPartANOVA[[i]]<-as.data.frame(varPartANOVA[[i]])
    varPartANOVA[[i]]$SumCellComp<-rowSums(varPartANOVA[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
	
	varPartIDOL[[i]]<-as.data.frame(varPartANOVA[[i]])
    varPartIDOL[[i]]$SumCellComp<-rowSums(varPartANOVA[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
	
	varPartANOVA[[i]]<-as.data.frame(varPartANOVA[[i]])
    varPartANOVA[[i]]$SumCellComp<-rowSums(varPartANOVA[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
}


par(mfrow = c(2,4))

	plotVarPart( varPartANOVA[[i]] )
}

plotVarPart( pcaPartANOVA[[1]] )

