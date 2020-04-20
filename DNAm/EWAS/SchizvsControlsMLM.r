## EWAS of schizophrenia
## uses a multi level model with random effect for individual & brain bank
## runs in parallel

runEWAS<-function(row,pheno){
	pheno<-cbind(row,pheno)
	model<-lmer(row ~ Phenotype + CCDNAmAge + Sex + Cell.type + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	ctNull<-lmer(row ~ Phenotype + CCDNAmAge + Sex  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	modelInt<-lmer(row ~ Phenotype + CCDNAmAge + Sex + Cell.type + Phenotype * Cell.type + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	outtab<-c(summary(model)$coefficients["PhenotypeControl",c(1,2,5)],
	anova(model, ctNull)[2,8],
	summary(model)$coefficients["Cell.typeNeuN +ve",c(1,2,5)],
	summary(model)$coefficients["Cell.typeSox10 +ve",c(1,2,5)],
	unlist(summary(model)$varcor),
	anova(modelInt,model)[2,8],
	summary(modelInt)$coefficients["PhenotypeControl",c(1,2,5)],
	summary(modelInt)$coefficients["Cell.typeNeuN +ve",c(1,2,5)],
	summary(modelInt)$coefficients["Cell.typeSox10 +ve",c(1,2,5)],
	summary(modelInt)$coefficients["PhenotypeControl:Cell.typeNeuN +ve",c(1,2,5)],
	summary(modelInt)$coefficients["PhenotypeControl:Cell.typeSox10 +ve",c(1,2,5)])
	return(outtab)
}


library(lme4)
library(lmerTest)
library(doParallel)

args<-commandArgs(trailingOnly = TRUE)
source(args[1])

setwd(dataDir)

load(normData)

## remove total samples and cell types with less than 10 samples
pheno<-pheno[which(pheno$Cell.type != "Total"),]
nSample<-table(pheno$Cell.type)
pheno<-pheno[pheno$Cell.type %in% names(nSample[which(nSample > 9)]),]

celltypenormbeta<-celltypenormbeta[,pheno$Basename]

## create parallel environment

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "lmer"))

#parRapply(cl, celltypenormbeta[1:10,], runEWAS, pheno)

outtab<-foreach(i=1:nrow(celltypenormbeta), .combine = "cbind") %dopar% runEWAS(celltypenormbeta[i,], pheno)
outtab<-t(outtab)
rownames(outtab)<-rownames(celltypenormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_P", "NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P", "RE_ID_SD", "RE_BrainBank_SD", "CellType_SCZ_P", "SCZ_coeff", "SCZ_SE", "SCZ_P", "NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P")

save(outtab, file = "Analysis/Schizophrenia/MLM.rdata")

stopCluster(cl)