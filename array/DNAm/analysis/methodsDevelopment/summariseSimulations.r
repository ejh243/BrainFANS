##---------------------------------------------------------------------#
##
## Title: Summarise Cell-specific EWAS Simulation results
##
## Purpose of script: collate and summarise result from EWAS simulations
## to compare different models
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

pThres<-9e-8
pThres.opts<-c(9e-8, 1e-7, 1e-6, 1e-5)

library(paletteer)
methodCols <- paletteer_d("ggsci::category10_d3")[1:4]
names(methodCols)<-c("ctLR", "allLR", "MER", "CRR")
methodLabels<-c("allLR: Main Effect", "allLR: Interaction", 
    "MER: Main Effect", "MER: Interaction", 
    "CRR: Main Effect", "CRR: Interaction", 
    "ctLR: DoubleNeg", "ctLR: NeunPos", "ctLR: Sox10Pos")
names(methodLabels)<-c("LM_ME", "LM_Int", "MLM_ME", "MLM_Int", "CRR_ME", "CRR_Int", "LM_Double-", "LM_NeuN+", "LM_Sox10+")


#----------------------------------------------------------------------#
# DEFINE FUNCTIONS
#----------------------------------------------------------------------#

countNSig<-function(data, thres){
	return(length(which(data < thres)))
}

countNSigMatrix<-function(matrix, thres){
	return(apply(matrix, 2, countNSig, thres))
}

addListNames<-function(orgList){
	newList <-lapply(names(orgList), function(name) {
		mat <- orgList[[name]]
		return(cbind(mat, pThres = name))
		})
	return(newList)
}

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape)


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)

## Simulation Scenario 1: across cell types simulate null ewas compare three analytical models

outFiles<-list.files(pattern="nullSimulations_Chunk")
outFiles<-outFiles[!duplicated(sub("(_MeanDiff.*)", "", outFiles))]
allSimComb<-NULL
for(each in outFiles){
	load(each)
	allSimComb<-cbind(allSimComb, nullSim)
}

#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#

fig0a<-list()
fig0b<-list()
xMax<--log10(min(allSimComb))
# look at p-value distribution to determine multiple testing burden
index<-1
quantileFullSimMinP<-rep(NA, 9)
for(each in c("LM_ME", "LM_Int", "MLM_ME", "MLM_Int", "CRR_ME", "CRR_Int", "LM_Double-", "LM_NeuN+", "LM_Sox10+")){
	titleText<-methodLabels[each]
    colKeep <- which(colnames(allSimComb) == each)
	## plot p-value distribution
	fig0a[[index]]<-melt(-log10(allSimComb[,each])) %>%
	ggplot(aes(x=value)) + 
	 ylab("nSimulations") + 
	 geom_histogram(bins = 50, boundary = 0, fill=methodCols[strsplit(titleText, ":")[[1]][1]]) + 
     xlab("-log10(P)") + labs(title = titleText) +
     xlim(0, xMax)

	# Finding the threshold (5th percentile of minimum p-values)
	FullSimMinP<-apply(allSimComb[,colKeep], 2, min)
	quantileFullSimMinP[index]<-quantile(FullSimMinP, 0.05)
	
	fig0b[[index]]<-ggplot(-log10(data.frame(FullSimMinP)), aes(x=FullSimMinP)) + 
	 geom_histogram(bins = 20, boundary = 0,  fill=methodCols[strsplit(titleText, ":")[[1]][1]]) + 
	 xlab("Maximum -log10(P)") + 
	 ylab("nSimulations")+ 
	 geom_vline(xintercept = -log10(quantileFullSimMinP[index])) + labs(title = titleText)+
     xlim(0, xMax)

	index<-index+1
}
fig0a<-ggarrange(plotlist=fig0a[c(seq(1,5,2), seq(2,6,2), 7:9)], nrow = 3, ncol = 3)
fig0b<-ggarrange(plotlist=fig0b[c(seq(1,5,2), seq(2,6,2), 7:9)], nrow = 3, ncol = 3)

ggsave(file.path(dataDir, "Plots", "HistogramNullEWASPvalue.pdf"), plot = fig0a, height = 8, width = 12) 
ggsave(file.path(dataDir, "Plots", "HistogramNullEWASMinP.pdf"), plot = fig0b, height = 8, width = 12) 

pos <- position_dodge(0.9)
## count number of significant associations per simulation
for(pThres in pThres.opts){
	nDMPs<-countNSigMatrix(allSimComb, pThres)
	dmpSumStats<-data.frame(matrix(unlist(strsplit(names(nDMPs), "_")), ncol = 2, byrow = TRUE), nDMPs)
	colnames(dmpSumStats)[c(1:2)]<-c("Method", "DMPtype")

	dmpSumStats$Method[dmpSumStats$DMPtype %in% c("Double-", "NeuN+", "Sox10+")]<-"ctLR"
	dmpSumStats$Method<- gsub("MLM", "MER", dmpSumStats$Method)
	dmpSumStats$Method<- gsub("LM", "allLR", dmpSumStats$Method)
	dmpSumStats$Method<-factor(dmpSumStats$Method, levels = c("ctLR", "allLR", "MER", "CRR"))
	dmpSumStats$DMPtype<-factor(dmpSumStats$DMPtype, levels = c("Double-", "NeuN+", "Sox10+", "ME", "Int"))
	#boxplot(nDMPs ~ as.factor(names(nDMPs)), ylab = "Number of DMPs", xlab = "Analytical model")

	y_lim<-range(nDMPs)
	fig1a <- ggplot(dmpSumStats, aes(x=Method, y=nDMPs, fill = Method))  +
	geom_violin(position = pos, scale = 'width')+
	stat_summary(fun = "mean", 
				geom = "point", 
				position = pos) + 
	scale_fill_manual(values = methodCols, labels = methodLabels, name = "Method") +
	ylim(y_lim)  +   
	theme(text = element_text(size = 20))  + 
	ylab("Number of false positives")    + 
	xlab("Method")  +
	facet_wrap(vars(DMPtype), nrow = 1, 
	labeller = labeller(DMPtype = c("ME" = "Main effect", "Int" = "Interaction", "Double-" = "DoubleNeg", "NeuN+" = "NeuNPos", "Sox10+" = "Sox10Pos")), scales = "free_x") +
	scale_y_continuous(trans='log10')

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotnDMPsNullEWASPThres", pThres, ".pdf")), plot = fig1a, height = 4, width = 12) 

	write.csv(aggregate(nDMPs ~ Method, dmpSumStats, FUN = function(x) c(mean = mean(x), sd = sd(x))),
	    file.path(dataDir, "Tables", paste0("SummaryStatsnDMPsByMethodNullEWASPThres", pThres, ".csv")))
	write.csv(aggregate(nDMPs ~ Method*DMPtype, dmpSumStats, FUN = function(x) c(mean = mean(x), sd = sd(x))),
	    file.path(dataDir, "Tables", paste0("SummaryStatsnDMPsByMethodxTermNullEWASPThres", pThres, ".csv")))
}

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

## Simulation Scenario 2: introduce true effects, prespecified to be either cell-type specific or common to all cell-types

outFiles<-list.files(pattern="nSigSimulateTrueEffects_Chunk")
outFiles<-outFiles[grep("pThres", outFiles)]
sumSimComb<-lapply(pThres.opts, function(x) NULL)
for(each in outFiles){
	meanDiff<-unlist(strsplit(strsplit(each, "_MeanDiff")[[1]][2], "\\_pThres_1e-05.rdata"))
	load(each)
	sumSim<-lapply(sumSim, function(mat) {
        mat<-cbind(as.numeric(meanDiff), mat)
        colnames(mat)[1]<-"MeanDiff"
        return(mat)
    })
	
	sumSimComb<-mapply(rbind, sumSimComb, sumSim, SIMPLIFY = FALSE)
}

names(sumSimComb) <- names(sumSim)

sumSimComb<-lapply(sumSimComb, as.data.frame)
sumSimComb<-lapply(sumSimComb, function(mat) {
        mat$propCS <- mat$nCTspecific/mat$nProbes
        return(mat)
    })

## compare performance at detecting true associations
## ME analysis should pick up all cell consistent effects and potentially cell specific effects also
## Int analysis should pick up cell-specific effects

sumSimComb <- lapply(sumSimComb, replace_na, list(0))

# fix NAs in simulations with no cell specific effects

sumSimComb <- lapply(sumSimComb, function(mat){
	mat[which(mat$nCTspecific == 0), "Double-_nCSTrueDMPs"] <- 0
	mat[which(mat$nCTspecific == 0),"NeuN+_nCSTrueDMPs"] <- 0
	mat[which(mat$nCTspecific == 0),"Sox10+_nCSTrueDMPs"] <- 0
	return(mat)
})

## calculate True positive rate for all DMPs (regardless of type) and regardless of which term in the model
tpRate<- lapply(sumSimComb, function(mat){
	tpRate<-cbind(mat$MeanDiff, (mat$LM_ME_nSigTrueDMPs + mat$LM_Int_nSigTrueDMPs)/(mat$nProbes), 
		(mat$MLM_ME_nSigTrueDMPs + mat$MLM_Int_nSigTrueDMPs)/(mat$nProbes),
		(mat$CRR_ME_nSigTrueDMPs + mat$CRR_Int_nSigTrueDMPs)/(mat$nProbes),
		mat$"LM_Double-_nSigTrueDMPs"/(mat$nProbes - mat$nCTspecific + mat$"Double-_nCSTrueDMPs"), 
		mat$"LM_NeuN+_nSigTrueDMPs"/(mat$nProbes - mat$nCTspecific + mat$"NeuN+_nCSTrueDMPs"), 
		mat$"LM_Sox10+_nSigTrueDMPs"/(mat$nProbes - mat$nCTspecific + mat$"Sox10+_nCSTrueDMPs"))
	colnames(tpRate)<-c("MeanDiff", "allLR_allLR", "MER_MER", "CRR_CRR", "ctLR_Double-", "ctLR_NeuN+", "ctLR_Sox10+")
	return(as.data.frame(tpRate))
})

tpRate <- addListNames(tpRate)
tpRate<-do.call("rbind", tpRate)


## calculate True positive rate for consistent effects
tpRate.CCon<-lapply(sumSimComb, function(mat){
	tpRate.CCon<-cbind(mat$MeanDiff, (mat$LM_ME_nSigTrueCommon)/(mat$nProbes-mat$nCTspecific), 
		(mat$LM_Int_nSigTrueCommon)/(mat$nProbes-mat$nCTspecific), 
		(mat$MLM_ME_nSigTrueCommon)/(mat$nProbes-mat$nCTspecific),
		(mat$MLM_Int_nSigTrueCommon)/(mat$nProbes-mat$nCTspecific), 
		(mat$CRR_ME_nSigTrueCommon)/(mat$nProbes-mat$nCTspecific),
		(mat$CRR_Int_nSigTrueCommon)/(mat$nProbes-mat$nCTspecific),
		(mat$"LM_Double-_nSigTrueCommon")/(mat$nProbes-mat$nCTspecific),		
		(mat$"LM_NeuN+_nSigTrueCommon")/(mat$nProbes-mat$nCTspecific),		
		(mat$"LM_Sox10+_nSigTrueCommon")/(mat$nProbes-mat$nCTspecific)		
		)
	colnames(tpRate.CCon)<-c("MeanDiff", "allLR_ME","allLR_Int", 
	"MER_ME", "MER_Int", 
	"CRR_ME", "CRR_Int", 
	"ctLR_Double-", "ctLR_NeuN+", "ctLR_Sox10+")
	return(as.data.frame(tpRate.CCon))
})
tpRate.CCon <- addListNames(tpRate.CCon)
tpRate.CCon<-do.call("rbind", tpRate.CCon)


## calculate True positive rate for ct specific effects
tpRate.CT<-lapply(sumSimComb, function(mat){
	tpRate.CT<-cbind(mat$MeanDiff, (mat$LM_ME_nSigTrueCS)/(mat$nCTspecific),
		(mat$LM_Int_nSigTrueCS)/(mat$nCTspecific), 
		(mat$MLM_ME_nSigTrueCS)/(mat$nCTspecific),
		(mat$MLM_Int_nSigTrueCS)/(mat$nCTspecific),
		(mat$CRR_ME_nSigTrueCS)/(mat$nCTspecific),
		(mat$CRR_Int_nSigTrueCS)/(mat$nCTspecific),
		(mat$"Double-_nSigTrueCS")/mat$"Double-_nCSTrueDMPs", 
		(mat$"NeuN+_nSigTrueCS")/mat$"NeuN+_nCSTrueDMPs",
		(mat$"Sox10+_nSigTrueCS")/mat$"Sox10+_nCSTrueDMPs"
		)
	colnames(tpRate.CT)<-c("MeanDiff", "allLR_ME","allLR_Int", 
	"MER_ME", "MER_Int", 
	"CRR_ME", "CRR_Int", 
	"ctLR_Double-", "ctLR_NeuN+", "ctLR_Sox10+")
	return(as.data.frame(tpRate.CT))
})
tpRate.CT <- addListNames(tpRate.CT)
tpRate.CT<-do.call("rbind", tpRate.CT)

## calculate false positive rate
fpRate<-lapply(sumSimComb, function(mat){
	fpRate<-cbind(mat$MeanDiff, (mat$LM_ME_nSigOther)/(mat$LM_ME_TotSig),
		(mat$LM_Int_nSigOther)/(mat$LM_Int_TotSig), 
		(mat$MLM_ME_nSigOther)/(mat$MLM_ME_TotSig),
		(mat$MLM_Int_nSigOther)/(mat$MLM_Int_TotSig),
		(mat$CRR_ME_nSigOther)/(mat$CRR_ME_TotSig),
		(mat$CRR_Int_nSigOther)/(mat$CRR_Int_TotSig),
		(mat$"LM_Double-_nSigOther")/(mat$"LM_Double-_TotSig"),		
		(mat$"LM_NeuN+_nSigOther")/(mat$"LM_NeuN+_TotSig"),		
		(mat$"LM_Sox10+_nSigOther")/(mat$"LM_Sox10+_TotSig")
		)
	colnames(fpRate)<-c("MeanDiff", "allLR_ME","allLR_Int", 
	"MER_ME", "MER_Int", 
	"CRR_ME", "CRR_Int", 
	"ctLR_Double-", "ctLR_NeuN+", "ctLR_Sox10+")
	return(as.data.frame(fpRate))
})
fpRate <- addListNames(fpRate)
fpRate<-do.call("rbind", fpRate)


#----------------------------------------------------------------------#
# SUMMARY STATISTICS
#----------------------------------------------------------------------#

longDat<-as.data.frame(tpRate) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	melt(id = c("MeanDiff", "pThres")) %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "allLR","MER", "CRR"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR")))

write.csv(aggregate(value ~ ., longDat, FUN = function(x) c(mean = mean(x), sd = sd(x))),
file.path(dataDir, "Tables", "SummaryStatisticsCTEWASTPRatesAllDMPs.csv"))

longDat<-as.data.frame(tpRate.CCon) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0))%>% 
	melt(id = c("MeanDiff", "pThres")) %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR")))

write.csv(aggregate(value ~ ., longDat, FUN = function(x) c(mean = mean(x), sd = sd(x))),
file.path(dataDir, "Tables", "SummaryStatisticsCTEWASTPRatesConsistentDMPs.csv"))

longDat<-as.data.frame(tpRate.CT) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0))%>% 
	melt(id = c("MeanDiff", "pThres")) %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR")))

write.csv(aggregate(value ~ ., longDat, FUN = function(x) c(mean = mean(x), sd = sd(x))),
file.path(dataDir, "Tables", "SummaryStatisticsCTEWASTPRatesCellTypeSpecificDMPs.csv"))

longDat<-as.data.frame(fpRate) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0))%>% 
	melt(id = c("MeanDiff", "pThres")) %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR")))

write.csv(aggregate(value ~ ., longDat, FUN = function(x) c(mean = mean(x), sd = sd(x))),
file.path(dataDir, "Tables", "SummaryStatisticsCTEWASFPRates.csv"))

#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#

# Initially just consider standard EWAS significance

for(meanDiff in c(0.02,0.05)){
	fig2a<-list()
	fig2a[[1]] <- as.data.frame(tpRate) %>% 
	filter(MeanDiff == meanDiff & pThres == pThres) %>% 
	select(allLR_allLR:"ctLR_Sox10+") %>% 
	replace_na(list("ctLM_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	melt() %>% separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "allLR","MER", "CRR"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>%
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + 
	  scale_fill_manual(values = methodCols, labels = methodLabels, name = "Method")+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) + theme(text = element_text(size = 20))  + 
				   ylab("True positive rate")    + 
				   xlab("Method") + 
				   labs(title = "All DMPs")
				   
	fig2a[[2]] <- as.data.frame(tpRate.CCon) %>% 
	filter(MeanDiff == meanDiff & pThres == pThres) %>% 
	select(allLR_ME:"ctLR_Sox10+") %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0))%>% 
	melt() %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels, name = "Method")+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + ylab("True positive rate")    + xlab("Method") +labs(title = "Common DMPs")

	fig2a[[3]] <- as.data.frame(tpRate.CT) %>% 
	filter(MeanDiff == meanDiff & pThres == pThres) %>% 
	select(allLR_ME:"ctLR_Sox10+") %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	melt() %>% separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels)+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + 
				   ylab("True positive rate")    + xlab("Method") +labs(title = "Cell-specific DMPs")


	fig2a<-ggarrange(plotlist=fig2a, nrow = 1, ncol = 3, common.legend = TRUE)

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotCTEWASTPFPRatesMeanDiff", meanDiff, ".pdf")), 
	plot = fig2a, height = 6, width = 18) 

}


# Look at effect of different p value thresholds

for(meanDiff in c(0.02,0.05)){
	fig3a <- as.data.frame(tpRate) %>% 
	filter(MeanDiff == meanDiff) %>% 
	select(allLR_allLR:pThres) %>% 
	replace_na(list("ctLM_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	melt(id = "pThres") %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "allLR","MER", "CRR"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>%
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + 
	  scale_fill_manual(values = methodCols, labels = methodLabels, name = "Method")+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) + theme(text = element_text(size = 20))  + 
	  ylab("True positive rate")    + 
	  xlab("Method") + 
	  labs(title = "All DMPs") +
	  facet_wrap(vars(pThres))

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotCTEWASTPRateAllDMPsMeanDiff", meanDiff, "FacetPThreshold.pdf")), 
	plot = fig3a, height = 12, width = 12) 
				   
	fig3b<- as.data.frame(tpRate.CCon) %>% 
	filter(MeanDiff == meanDiff) %>% 
	select(allLR_ME:pThres) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0))%>% 
	melt(id = "pThres") %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels, name = "Method")+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + 
		ylab("True positive rate")    + 
		xlab("Method") +
		labs(title = "Common DMPs") +
	  facet_wrap(vars(pThres))

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotCTEWASTPRateConsistentDMPsMeanDiff", meanDiff, "FacetPThreshold.pdf")), 
	plot = fig3b, height = 12, width = 12) 

	fig3c <- as.data.frame(tpRate.CT) %>% 
	filter(MeanDiff == meanDiff) %>% 
	select(allLR_ME:pThres) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	melt(id = "pThres") %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels)+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + 
		ylab("True positive rate")    + 
		xlab("Method") +
		labs(title = "Cell-specific DMPs") +
	  facet_wrap(vars(pThres))

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotCTEWASTPRateCellSpecificDMPsMeanDiff", meanDiff, "FacetPThreshold.pdf")), 
	plot = fig3c, height = 12, width = 12) 
				   
	fig3d <- as.data.frame(fpRate)%>% 
	filter(MeanDiff == meanDiff) %>% 
	select(allLR_ME:pThres) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0))%>% 
	melt(id = "pThres") %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% 
	mutate(Regression = factor(Regression, levels = c("ctLR", "allLR", "MER", "CRR"))) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels)+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + 
		ylab("False positive rate")    + 
		xlab("Method") +
		labs(title = "") +
	  facet_wrap(vars(pThres))

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotCTEWASFPRateMeanDiff", meanDiff, "FacetPThreshold.pdf")), 
	plot = fig3d, height = 12, width = 12) 
}



fig3a<-list()
tpRate<-as.data.frame(tpRate)
tpRate$nProbes<-unlist(lapply(sumSimComb, "[", "nProbes"))
tpRate$propCS<-unlist(lapply(sumSimComb, "[", "propCS"))


for(meanDiff in c(0.02,0.05)){
	fig3a<-list()
	fig3a[[1]] <- as.data.frame(tpRate) %>% 
	filter(MeanDiff == meanDiff & pThres == pThres) %>% 
	select(allLR_allLR:nProbes) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	group_by(nProbes) %>% 
	summarise(ctLR_DoubleNeg = mean(`ctLR_Double-`), ctLR_NeuNPos = mean(`ctLR_NeuN+`), ctLR_Sox10Pos = mean(`ctLR_Sox10+`),
	    allLR_allLR = mean(allLR_allLR), MER_MER = mean(MER_MER), CRR_CRR = mean(CRR_CRR)) %>%
  pivot_longer(
    cols = -nProbes,
    names_to = "Measurement",
    values_to = "Value") %>% 
	ggplot( aes(x=nProbes, y=Value, colour = Measurement))  +
	  geom_line(size = 2) + 
				   ylab("True positive rate")    + 
				   xlab("Number of true DMPs") + 
				   labs(title = "All DMPs")


## as function of proportion of cell-specific effects
fig3a[[2]]<-as.data.frame(tpRate) %>% 
	filter(MeanDiff == meanDiff & pThres == pThres) %>% 
	select(allLR_allLR:"ctLR_Sox10+", propCS) %>% 
	replace_na(list("ctLR_Double-" = 0, "ctLR_NeuN+" = 0, "ctLR_Sox10+" = 0)) %>% 
	group_by(propCS) %>% 
	summarise(ctLR_DoubleNeg = mean(`ctLR_Double-`), ctLR_NeuNPos = mean(`ctLR_NeuN+`), ctLR_Sox10Pos = mean(`ctLR_Sox10+`),
	    allLR_allLR = mean(allLR_allLR), MER_MER = mean(MER_MER), CRR_CRR = mean(CRR_CRR)) %>%
  pivot_longer(
    cols = -propCS,
    names_to = "Measurement",
    values_to = "Value") %>% 
	ggplot( aes(x=propCS, y=Value, colour = Measurement))  +
	  geom_line(size = 2) + 
				   ylab("True positive rate")    + 
				   xlab("Proportion DMPs cell-type specific") + 
				   labs(title = "All DMPs")



	fig3a<-ggarrange(plotlist=fig3a, nrow = 1, ncol = 2, common.legend = TRUE)
	ggsave(file.path(dataDir, "Plots", paste0("LineGraphCTEWASTPFPRatesMeanDiff", meanDiff, ".pdf")), plot = fig3a, height = 6, width = 12) 
}

