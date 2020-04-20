## EWAS of schizophrenia
## uses a standard linear regression model 

args<-commandArgs(trailingOnly = TRUE)
source(args[1])

setwd(dataDir)

load(normData)

## remove total samples and cell types with less than 10 samples
pheno<-pheno[which(pheno$Cell.type != "Total"),]
nSample<-table(pheno$Cell.type)
pheno<-pheno[pheno$Cell.type %in% names(nSample[which(nSample > 9)]),]

celltypenormbeta<-celltypenormbeta[,pheno$Basename]

outtab<-matrix(data = NA, ncol = 26, nrow = nrow(celltypenormbeta))
rownames(outtab)<-rownames(celltypenormbeta)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", "CellType_P", "NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P", "CellType_SCZ_P", "SCZ_coeff", "SCZ_SE", "SCZ_P", "NeuN_coeff", "NeuN_SE", "NeuN_P", "SOX10_coeff", "SOX10_SE", "SOX10_P", "NeuN_SCZ_coeff", "NeuN_SCZ_SE","NeuN_SCZ_P", "SOX10_SCZ_coeff", "SOX10_SCZ_SE","SOX10_SCZ_P")

for(i in 1:nrow(celltypenormbeta)){

	model<-lm(celltypenormbeta[i,] ~ pheno$Phenotype + pheno$CCDNAmAge + pheno$Sex + pheno$Cell.type + pheno$Tissue.Centre)
	ctNull<-lm(celltypenormbeta[i,] ~ pheno$Phenotype + pheno$CCDNAmAge + pheno$Sex  + pheno$Tissue.Centre)
	modelInt<-lm(celltypenormbeta[i,] ~ pheno$Phenotype + pheno$CCDNAmAge + pheno$Sex + pheno$Cell.type + pheno$Phenotype * pheno$Cell.type + pheno$Tissue.Centre)
	outtab[i,1:3]<-summary(model)$coefficients["pheno$PhenotypeControl",c(1,2,4)]
	outtab[i,4]<-anova(model, ctNull)[2,6]
	outtab[i,5:7]<-summary(model)$coefficients["pheno$Cell.typeNeuN +ve",c(1,2,4)]
	outtab[i,8:10]<-summary(model)$coefficients["pheno$Cell.typeSox10 +ve",c(1,2,4)]
	outtab[i,11]<-anova(modelInt,model)[2,6]
	outtab[i,12:14]<-summary(modelInt)$coefficients["pheno$PhenotypeControl",c(1,2,4)]
	outtab[i,15:17]<-summary(modelInt)$coefficients["pheno$Cell.typeNeuN +ve",c(1,2,4)]
	outtab[i,18:20]<-summary(modelInt)$coefficients["pheno$Cell.typeSox10 +ve",c(1,2,4)]
	outtab[i,21:23]<-summary(modelInt)$coefficients["pheno$PhenotypeControl:pheno$Cell.typeNeuN +ve",c(1,2,4)]
	outtab[i,24:26]<-summary(modelInt)$coefficients["pheno$PhenotypeControl:pheno$Cell.typeSox10 +ve",c(1,2,4)]	
}

save(outtab, file = "Analysis/Schizophrenia/LM.rdata")