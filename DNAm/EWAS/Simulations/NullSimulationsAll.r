## compare EWAS methods on simulated data
## randomly select individuals as cases and controls to match distribution per brain bank
## compare 3 models analyses data using data all cell types

runEWAS<-function(row,pheno, status){

	modelLM<-lm(row ~ status * pheno$Cell.type + pheno$CCDNAmAge + pheno$Sex + pheno$Tissue.Centre)
	nullLM<-lm(row ~ status + pheno$Cell.type + pheno$CCDNAmAge + pheno$Sex + pheno$Tissue.Centre)

	data<-cbind(row,status, pheno)
	modelMLM<-lmer(row ~ status * Cell.type + CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)
	nullMLM<-lmer(row ~ status + Cell.type + CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = pheno)

	p.df <- pdata.frame(data.frame("meth" = row, "phenotype" = status, "age" = pheno$CCDNAmAge, "sex" = pheno$Sex, "cell.type" = pheno$Cell.type, "brain.bank" = pheno$Tissue.Centre, "id" = pheno$Indidivual.ID), index = c("id", "cell.type"), drop.index = F)
		
	modelp <- plm(meth ~ phenotype * cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	nullCRR <- plm(meth ~ phenotype + cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	# compute Stata like df-adjustment
	G <- length(unique(p.df$id))
	N <- length(p.df$id)
	dfa <- (G/(G - 1)) * (N - 1)/modelp$df.residual

	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(modelp, type = "HC0", cluster = "group", adjust = T)
	modelCRR<-coeftest(modelp, vcov = firm_c_vcov)


	return(c(summary(modelLM)$coefficients["status1", 4], anova(modelLM, nullLM)[2,6], 
	summary(modelMLM)$coefficients["status1", 5], anova(modelMLM, nullMLM)[2,8],
	 modelCRR["phenotype1",4],waldtest(modelp, nullCRR, vcov = firm_c_vcov, test = "F")[2,4]))
}

library(lme4)
library(lmerTest)
library(plm)
library(lmtest)
#library(GenABEL)
library(doParallel)

nSim<-10
nChunk<-10

args<-commandArgs(trailingOnly = TRUE)
source(args[1])

setwd(dataDir)

load(normData)

## remove total samples and cell types with less than 10 samples
pheno<-pheno[which(pheno$Cell.type != "Total"),]
nSample<-table(pheno$Cell.type)
pheno<-pheno[pheno$Cell.type %in% names(nSample[which(nSample > 9)]),]

celltypenormbeta<-celltypenormbeta[,pheno$Basename]


nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "lmer", "pdata.frame", "plm", "vcovHC", "coeftest", "waldtest"))


## first randomly assign cases control status, do we see inflation?
for(chunk in nChunk){
	allSim<-NULL
	for(simNum in 1:nSim){
		status<-rep(0,nrow(pheno))

		for(centre in unique(pheno$Tissue.Centre)){
			ids<-unique(pheno$Indidivual.ID[which(pheno$Tissue.Centre == centre)])
			cases<-sample(ids, floor(length(ids)/2))
			status[pheno$Indidivual.ID %in% cases]<-1
		}
		status<-as.factor(status)
		outtab<-foreach(i=1:nrow(celltypenormbeta), .combine = "cbind") %dopar% runEWAS(celltypenormbeta[i,], pheno, status)
		outtab<-t(outtab)

		rownames(outtab)<-rownames(celltypenormbeta)
		colnames(outtab)<-c(paste("LM", c("ME", "Int"), sep = "_"),paste("MLM", c("ME", "Int"), sep = "_"),paste("CRR", c("ME", "Int"), sep = "_"))
		allSim<-cbind(allSim, outtab)
	}

	save(allSim, file = paste0("nullSimulations_Chunk", chunk, ".rdata"))
}