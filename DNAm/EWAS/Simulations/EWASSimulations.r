## simulate EWAS with induced changes between cases and controls
## set proportion of differences cell-type specific vs universal
## compare 3 models
## run as chunks


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

estlambda<-function(pvals){
	z = qnorm(pvals / 2)
	## calculates lambda
	lambda = round(median(z^2) / 0.454, 3)
	return(lambda)
}

library(lme4)
library(lmerTest)
library(plm)
library(lmtest)
#library(GenABEL)
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

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "lmer", "pdata.frame", "plm", "vcovHC", "coeftest", "waldtest"))

cellTypes<-unique(pheno$Cell.type)

## use chunks to save intermitant results in case script does not complete
nSim<-10
nChunk<-10
nSig.options<-c(10,100,1000)
propCS.options<-seq(0,1,0.1)
sigEffect<-0.05
for(chunk in nChunk){
	sumSim<-matrix(data = NA, nrow = nSim*length(nSig.options)*length(propCS.options), ncol = 2+(6*6))
	colnames(sumSim)<-c("nProbes", "nCTspecific", paste("LM_ME", c("TotSig", "nTruePos", "nFalsePos","nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
	paste("LM_Int", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
	paste("MLM_ME", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
	paste("MLM_Int", c("TotSig", "nTruePos", "nFalsePos","nTruePosCS", "nFalsePosCS",  "lambda"), sep = "_"),
	paste("CRR_ME", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"),
	paste("CRR_Int", c("TotSig", "nTruePos", "nFalsePos", "nTruePosCS", "nFalsePosCS", "lambda"), sep = "_"))

	rowNum<-1
	for(simNum in 1:nSim){
		## randomly select samples to be cases, fix so numbers per brain bank match actual
		status<-rep(0,nrow(pheno))
		for(centre in unique(pheno$Tissue.Centre)){
			ids<-unique(pheno$Indidivual.ID[which(pheno$Tissue.Centre == centre)])
			cases<-sample(ids, floor(length(ids)/2))
			status[pheno$Indidivual.ID %in% cases]<-1
		}
		status<-as.factor(status)
		
		## First run the null EWAS for this simulation
		outtab.null<-foreach(i=1:nrow(celltypenormbeta), .combine = "rbind") %dopar% runEWAS(celltypenormbeta[i,], pheno, status)
		## Second only need to retest sites affects are induced at as non significiant probes are unaltered		
		for(nSig in nSig.options){
			for(propCS in propCS.options){
				sumSim[rowNum,1]<-nSig
				sumSim[rowNum,2]<-floor(nSig*propCS)
				## randomly select significant probes
				sigProbes<-sample(1:nrow(celltypenormbeta), nSig)

				## create matrix with effects to add to sample profiles
				diffs<-rnorm(nrow(pheno)*nSig, sigEffect, 0.005) ## sample effects for each sample
				diffs<-matrix(data = diffs, nrow = nSig, byrow = TRUE)
				diffs[,which(status == 0)]<-0
				
				## make some of these effects cell type specific
				if(propCS > 0){
					ctSpecific<-sample(1:nSig, floor(nSig*propCS))
					for(each in ctSpecific){
						# randomly select cell type to be significant in
						selectCT<-sample(cellTypes, 1)
						## set other cell types to have no effects
						diffs[each, !pheno$Cell.type %in% selectCT]<-0
					}
					ctSpecific<-sigProbes[ctSpecific]
				} else {
					ctSpecific<-NULL
				}				

				## create matrix of betas with induced effects
				testbetas<-celltypenormbeta[sigProbes,]+diffs

				outtab.sig<-foreach(i=1:nrow(testbetas), .combine = "rbind") %dopar% runEWAS(testbetas[i,], pheno, status)
				
				## merge signif results with null results to generate summary statistics
				outtab.sim<-outtab.null
				outtab.sim[sigProbes,]<-outtab.sig
				
				sumSim[rowNum,2+seq(1,6*6, 6)]<-colSums(outtab.sim < 9e-8)
				sumSim[rowNum,3+seq(1,6*6, 6)]<-colSums(outtab.sim[sigProbes,] < 9e-8)
				sumSim[rowNum,4+seq(1,6*6, 6)]<-colSums(outtab.sim[-sigProbes,] < 9e-8)
				## handle quirk of R converting 1 row matrix to vector!!
				if(length(ctSpecific) > 2){
					sumSim[rowNum,5+seq(1,6*6, 6)]<-colSums(outtab.sim[ctSpecific,] < 9e-8)
					sumSim[rowNum,6+seq(1,6*6, 6)]<-colSums(outtab.sim[-ctSpecific,] < 9e-8)
				} else {
					if(length(ctSpecific) == 1){
						sumSim[rowNum,5+seq(1,6*6, 6)]<-as.numeric(outtab.sim[ctSpecific,] < 9e-8)
						sumSim[rowNum,6+seq(1,6*6, 6)]<-colSums(outtab.sim[-ctSpecific,] < 9e-8)
					} else{
						sumSim[rowNum,5+seq(1,6*6, 6)]<-NA
						sumSim[rowNum,6+seq(1,6*6, 6)]<-NA
					}
				}
							
				sumSim[rowNum,7+seq(1,6*6, 6)]<-apply(outtab.sim, 2, estlambda)
				rowNum<-rowNum+1
			}
		}
	}

	save(sumSim, file = paste0("nSigSimulateTrueEffects_Chunk", chunk, ".rdata"))
}