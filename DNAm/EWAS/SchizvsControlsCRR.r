## EWAS of schizophrenia
## uses clustered robust regression

library(plm)
library(lmtest)

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

	p.df <- pdata.frame(data.frame("meth" = celltypenormbeta[i,], "phenotype" = pheno$Phenotype, "age" = pheno$CCDNAmAge, "sex" = pheno$Sex, "cell.type" = pheno$Cell.type, "brain.bank" = pheno$Tissue.Centre, "id" = pheno$Indidivual.ID), index = c("id", "cell.type"), drop.index = F)
	
	pmodel <- plm(meth ~ phenotype + age + sex + cell.type + brain.bank, data = p.df, model = "pooling")
	# compute Stata like df-adjustment
	G <- length(unique(p.df$id))
	N <- length(p.df$id)
	dfa <- (G/(G - 1)) * (N - 1)/pmodel$df.residual
	
	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(pmodel, type = "HC0", cluster = "group", adjust = T)
	model<-coeftest(pmodel, vcov = firm_c_vcov)

	ctNull<-plm(meth ~ phenotype + age + sex + brain.bank, data = p.df, model = "pooling")
	outtab[i,1:3]<-model["phenotypeControl",c(1,2,4)]
	outtab[i,4]<-waldtest(pmodel, ctNull, vcov = firm_c_vcov, test = "F")[2,4]
	outtab[i,5:7]<-model["cell.typeNeuN +ve",c(1,2,4)]
	outtab[i,8:10]<-model["cell.typeSox10 +ve",c(1,2,4)]	
	
	pmodelInt <- plm(meth ~ phenotype + age + sex + cell.type + brain.bank + phenotype*cell.type, data = p.df, model = "pooling")
	dfa <- (G/(G - 1)) * (N - 1)/pmodelInt$df.residual
	
	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(pmodelInt, type = "HC0", cluster = "group", adjust = T)
	modelInt<-coeftest(pmodelInt, vcov = firm_c_vcov)

	outtab[i,11]<-waldtest(pmodelInt,pmodel, vcov = firm_c_vcov, test = "F")[2,4]
	outtab[i,12:14]<-modelInt["phenotypeControl",c(1,2,4)]
	outtab[i,15:17]<-modelInt["cell.typeNeuN +ve",c(1,2,4)]
	outtab[i,18:20]<-modelInt["cell.typeSox10 +ve",c(1,2,4)]
	outtab[i,21:23]<-modelInt["phenotypeControl:cell.typeNeuN +ve",c(1,2,4)]
	outtab[i,24:26]<-modelInt["phenotypeControl:cell.typeSox10 +ve",c(1,2,4)]	
}

save(outtab, file = "Analysis/Schizophrenia/RobustClusteredErrors.rdata")