##---------------------------------------------------------------------#
##
## Title: Simulate Cell-specific EWAS
##
## Purpose of script: simulate EWAS with induced changes between cases and controls for both 
## cell-type specific & ubiquitous DMPs to compare analytical models
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runSimEWAS<-function(row,QCmetrics, status){

	# OLS across all CTs together
	modelLM<-lm(row ~ status * QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	nullLM<-lm(row ~ status + QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)

	#Mixed effects model
	data<-cbind(row,status, QCmetrics)
	modelMLM<-lmer(row ~ status * Cell.type + CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = QCmetrics)
	nullMLM<-lmer(row ~ status + Cell.type + CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = QCmetrics)

	#Clustered Robust Regression (CRR)
	p.df <- pdata.frame(data.frame("meth" = row, "phenotype" = status, "age" = QCmetrics$CCDNAmAge, "sex" = QCmetrics$Sex, "cell.type" = QCmetrics$Cell.type, "brain.bank" = QCmetrics$Tissue.Centre, "id" = QCmetrics$Indidivual.ID), index = c("id"), drop.index = F)
		
	modelp <- plm(meth ~ phenotype * cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	nullCRR <- plm(meth ~ phenotype + cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	# compute Stata like df-adjustment
	G <- length(unique(p.df$id))
	N <- length(p.df$id)
	dfa <- (G/(G - 1)) * (N - 1)/modelp$df.residual

	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(modelp, type = "HC0", cluster = "group", adjust = T)
	modelCRR<-coeftest(modelp, vcov = firm_c_vcov)

	# OLS within each cell type
	ctP<-NULL
	for(type in unique(QCmetrics$Cell.type)){
		modelCT<-lm(row ~ status + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre, subset = which(QCmetrics$Cell.type == type))
		ctP<-c(ctP, summary(modelCT)$coefficients["status1", 4])
	}

	#collate results for return
	return(c(summary(modelLM)$coefficients["status1", 4], anova(modelLM, nullLM)[2,6], 
	summary(modelMLM)$coefficients["status1", 5], anova(modelMLM, nullMLM)[2,8],
	 modelCRR["phenotype1",4],waldtest(modelp, nullCRR, vcov = firm_c_vcov, test = "F")[2,4], ctP))
}

runCTEWAS<-function(row,QCmetrics){

	modelLM<-lm(row ~ QCmetrics$Cell.type + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
	nullLM<-lm(row ~ QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)

	data<-cbind(row,QCmetrics)
	modelMLM<-lmer(row ~ Cell.type + CCDNAmAge + Sex + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = QCmetrics)
	nullMLM<-lmer(row ~ CCDNAmAge + Sex +  + (1 | Tissue.Centre)  + (1 | Indidivual.ID), REML = FALSE, data = QCmetrics)

	p.df <- pdata.frame(data.frame("meth" = row, "age" = QCmetrics$CCDNAmAge, "sex" = QCmetrics$Sex, "cell.type" = QCmetrics$Cell.type, "brain.bank" = QCmetrics$Tissue.Centre, "id" = QCmetrics$Indidivual.ID), index = c("id"), drop.index = F)
		
	modelp <- plm(meth ~ cell.type + age + sex  + brain.bank, data = p.df, model = "pooling")
	nullCRR <- plm(meth ~ age + sex  + brain.bank, data = p.df, model = "pooling")
	# compute Stata like df-adjustment
	G <- length(unique(p.df$id))
	N <- length(p.df$id)
	dfa <- (G/(G - 1)) * (N - 1)/modelp$df.residual

	# display with cluster VCE and df-adjustment
	firm_c_vcov <- dfa * vcovHC(modelp, type = "HC0", cluster = "group", adjust = T)
	modelCRR<-coeftest(modelp, vcov = firm_c_vcov)

	return(c(anova(modelLM, nullLM)[2,6], 
	anova(modelMLM, nullMLM)[2,8],
	waldtest(modelp, nullCRR, vcov = firm_c_vcov, test = "F")[2,4], 
	leveneTest(row ~ QCmetrics$Cell.type)[1,3],
	aggregate(row ~ QCmetrics$Cell.type, FUN = mean)$row,
    aggregate(row ~ QCmetrics$Cell.type, FUN = sd)$row))

}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(lme4)
library(lmerTest)
library(plm)
library(lmtest)
library(doParallel)
library(car)
library(cdegUtilities)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)

nSim<-10
nSig.options<-c(10,100,1000)
propCS.options<-seq(0,1,0.2)
pThres.opts<-c(9e-8, 1e-7, 1e-6, 1e-5)

dataDir <- args[1]
nChunk<-args[2]
sigEffect<-as.numeric(args[3])

set.seed(nChunk)

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis/methodsDevelopment")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 100 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 99)]),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

cellTypes<-unique(QCmetrics$Cell.type)

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores/2)
registerDoParallel(cl)
clusterExport(cl, list("runSimEWAS", "runCTEWAS", "lmer", "pdata.frame", "plm", "vcovHC", "coeftest", "waldtest", "leveneTest", "pThres.opts"))

#----------------------------------------------------------------------#
# RUN CELL TYPE EWAS
#----------------------------------------------------------------------#

# Run the ct EWAS for all probes 
ctEWAS<-foreach(i=1:nrow(celltypeNormbeta), .combine = "rbind") %dopar% runCTEWAS(celltypeNormbeta[i,], QCmetrics)
colnames(ctEWAS)<-c("LM_ctP", "MLM_ctP", "CRR_ctP", "LevenesP", outer(unique(QCmetrics$Cell.type), c("mean", "sd"), FUN = "paste", sep = "_"))

#----------------------------------------------------------------------#
# RUN SIMULATIONS
#----------------------------------------------------------------------#

sumSim<-matrix(data = NA, nrow = nSim*length(nSig.options)*length(propCS.options), ncol = 2+(6*9)+(2*3))
colnames(sumSim)<-c("nProbes", "nCTspecific", paste("LM_ME", c("TotSig", "nSigTrueDMPs", "nSigOther","nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"),
paste("LM_Int", c("TotSig", "nSigTrueDMPs", "nSigOther", "nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"),
paste("MLM_ME", c("TotSig", "nSigTrueDMPs", "nSigOther", "nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"),
paste("MLM_Int", c("TotSig", "nSigTrueDMPs", "nSigOther","nSigTrueCS", "nSigTrueCommon",  "lambda"), sep = "_"),
paste("CRR_ME", c("TotSig", "nSigTrueDMPs", "nSigOther", "nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"),
paste("CRR_Int", c("TotSig", "nSigTrueDMPs", "nSigOther", "nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"),
paste("LM", cellTypes[1], c("TotSig", "nSigTrueDMPs", "nSigOther","nSigTrueCS", "nSigTrueCommon",  "lambda"), sep = "_"),
paste("LM", cellTypes[2], c("TotSig", "nSigTrueDMPs", "nSigOther", "nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"),
paste("LM", cellTypes[3], c("TotSig", "nSigTrueDMPs", "nSigOther", "nSigTrueCS", "nSigTrueCommon", "lambda"), sep = "_"), 
outer(cellTypes, c("nCSTrueDMPs", "nSigTrueCS"), FUN = "paste", sep = "_"))

sumSim<-lapply(pThres.opts, function(x) sumSim)
names(sumSim)<-pThres.opts

rowNum<-1
nullSim<-NULL

commonDMPs<-NULL
ctDMPs<-NULL
fpDMPs<-NULL
for(simNum in 1:nSim){
	message("Simulation: ", simNum, " of ", nSim)

	# randomly select samples to be cases, fix so numbers per brain bank match actual
	status<-rep(0,nrow(QCmetrics))
	for(centre in unique(QCmetrics$Tissue.Centre)){
		ids<-unique(QCmetrics$Indidivual.ID[which(QCmetrics$Tissue.Centre == centre)])
		cases<-sample(ids, floor(length(ids)/2))
		status[QCmetrics$Indidivual.ID %in% cases]<-1
	}
	status<-as.factor(status)
	
	# Run the null EWAS for this simulation
	outtab.null<-foreach(i=1:nrow(celltypeNormbeta), .combine = "rbind") %dopar% runSimEWAS(celltypeNormbeta[i,], QCmetrics, status)
	
	rownames(outtab.null)<-rownames(celltypeNormbeta)
	colnames(outtab.null)<-c(paste("LM", c("ME", "Int"), sep = "_"),paste("MLM", c("ME", "Int"), sep = "_"),paste("CRR", c("ME", "Int"), sep = "_"), paste("LM", unique(QCmetrics$Cell.type), sep = "_"))
	nullSim<-cbind(nullSim, outtab.null)
	
	# Retest sites effects are induced at as non significiant probes are unaltered		
	for(nSig in nSig.options){
		for(propCS in propCS.options){
			for(pThres in pThres.opts){
				sumSim[[as.character(pThres)]][rowNum,1]<-nSig
				sumSim[[as.character(pThres)]][rowNum,2]<-floor(nSig*propCS)
			}
			# randomly select significant probes
			sigProbes<-sample(1:nrow(celltypeNormbeta), nSig)

			# create matrix with effects to add to sample profiles
			diffs<-rnorm(nrow(QCmetrics)*nSig, sigEffect, 0.005) 
			diffs<-matrix(data = diffs, nrow = nSig, byrow = TRUE)
			diffs[,which(status == 0)]<-0 # only add to cases
			
			# make soof these effects cell type specific
			if(propCS > 0){
				ctSpecific<-sample(1:nSig, floor(nSig*propCS))
				## record which ct specifc to
				affectedCT<-sapply(cellTypes,function(x) NULL)
				for(each in ctSpecific){
					# randomly select cell type to be significant in
					selectCT<-sample(cellTypes, 1)
					# set other cell types to have no effect
					diffs[each, !QCmetrics$Cell.type %in% selectCT]<-0
					# record 
					affectedCT[[selectCT]]<-c(affectedCT[[selectCT]], sigProbes[each])
				}
				ctSpecific<-sigProbes[ctSpecific]
			} else {
				ctSpecific<-NULL
			}				
			commonProbes<-sigProbes[!sigProbes %in% ctSpecific]
			
			# create matrix of betas with induced effects
			testbetas<-celltypeNormbeta[sigProbes,]+diffs

			outtab.sig<-foreach(i=1:nrow(testbetas), .combine = "rbind") %dopar% runSimEWAS(testbetas[i,], QCmetrics, status)
			
			# merge signif results with null results to generate summary statistics
			outtab.sim<-outtab.null
			outtab.sim[sigProbes,]<-outtab.sig
			
			for(pThres in pThres.opts){
				sumSim[[as.character(pThres)]][rowNum,2+seq(1,6*9,6)]<-colSums(outtab.sim < pThres)
				sumSim[[as.character(pThres)]][rowNum,3+seq(1,6*9,6)]<-colSums(outtab.sim[sigProbes,] < pThres)
				sumSim[[as.character(pThres)]][rowNum,4+seq(1,6*9,6)]<-colSums(outtab.sim[-sigProbes,] < pThres)
				
				# substract those significant with both main effect (ME) and interation term from ME count
				methodMaxP<-data.frame("LM" = apply(outtab.sim[,1:2], 1, max), "MLM" = apply(outtab.sim[,3:4], 1, max), "CRR" = apply(outtab.sim[,5:6], 1, max))
				sumSim[[as.character(pThres)]][rowNum,2+seq(1,6*6, 12)]<-sumSim[[as.character(pThres)]][rowNum,2+seq(1,6*6, 12)] - colSums(methodMaxP < pThres)
				sumSim[[as.character(pThres)]][rowNum,3+seq(1,6*6, 12)]<-sumSim[[as.character(pThres)]][rowNum,3+seq(1,6*6, 12)] - colSums(methodMaxP[sigProbes,] < pThres)
				sumSim[[as.character(pThres)]][rowNum,4+seq(1,6*6, 12)]<-sumSim[[as.character(pThres)]][rowNum,4+seq(1,6*6, 12)] - colSums(methodMaxP[-sigProbes,] < pThres)
				
				# handle quirk of R converting 1 row matrix to vector!!
				if(length(ctSpecific) > 1){
					sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*9,6)]<-colSums(outtab.sim[ctSpecific,] < pThres)	
					sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*6, 12)]<-sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*6, 12)] - colSums(methodMaxP[ctSpecific,] < pThres)
					
					# check ct specific DMPs detected in OLS within that cell type
					for(type in cellTypes){
						if(length(affectedCT[[type]]) > 0){
							sumSim[[as.character(pThres)]][rowNum,paste(type, c("nCSTrueDMPs", "nSigTrueCS"), sep = "_")]<-c(length(affectedCT[[type]]),sum(outtab.sim[affectedCT[[type]],which(colnames(outtab.sim) == paste("LM", type, sep = "_"))] < pThres))					 
						} else {
							sumSim[[as.character(pThres)]][rowNum,paste(type, c("nCSTrueDMPs", "nSigTrueCS"), sep = "_")]<-c(0,0)
						}
					}
					
				} else {
					if(length(ctSpecific) == 1){
						sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*9,6)]<-as.numeric(outtab.sim[ctSpecific,] < pThres)
						sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*6, 12)]<-sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*6, 12)] - as.numeric(methodMaxP[ctSpecific,] < pThres)
					} else{
						sumSim[[as.character(pThres)]][rowNum,5+seq(1,6*9,6)]<-0
					}
				}
				
				if(length(commonProbes) > 0){
					sumSim[[as.character(pThres)]][rowNum,6+seq(1,6*9,6)]<-colSums(outtab.sim[commonProbes,] < pThres)
					sumSim[[as.character(pThres)]][rowNum,6+seq(1,6*6, 12)]<-sumSim[[as.character(pThres)]][rowNum,6+seq(1,6*6, 12)] - colSums(methodMaxP[commonProbes,] < pThres)
				}else{
					sumSim[[as.character(pThres)]][rowNum,6+seq(1,6*9,6)]<-0
				}
							
				sumSim[[as.character(pThres)]][rowNum,7+seq(1,6*9,6)]<-apply(outtab.sim, 2, estlambda)
				
			}
			rowNum<-rowNum+1
			## collate sum stats of commonDMPs
			commonDMPs<-rbind(commonDMPs, cbind(outtab.sim[commonProbes,], ctEWAS[commonProbes,]))
			
			## collate sum stats of ctSpecificDMPs
			ctDMPs<-rbind(ctDMPs, cbind(outtab.sim[ctSpecific,], ctEWAS[ctSpecific,]))

			## collate sum stats of false positives
			fpIndex<-which(apply(outtab.sim, 1, min) < pThres)
			fpDMPs<-rbind(fpDMPs, cbind(outtab.sim[fpIndex,], ctEWAS[fpIndex,]))
		}
	}
}

save(nullSim, file = file.path(resPath, paste0("nullSimulations_Chunk", nChunk, "_MeanDiff",sigEffect, ".rdata")))
save(sumSim, file = file.path(resPath, paste0("nSigSimulateTrueEffects_Chunk", nChunk, "_MeanDiff",sigEffect, "_pThres_", pThres,  ".rdata")))
save(commonDMPs, ctDMPs, file = file.path(resPath, paste0("DMPsumStats_Chunk", nChunk, "_MeanDiff",sigEffect, "_pThres_", pThres,  ".rdata")))
