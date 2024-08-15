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

estlambda<-function(pvals){
	z = qnorm(pvals / 2)
	## calculates lambda
	lambda = round(median(z^2) / 0.454, 3)
	return(lambda)
}

estlambdaMatrix<-function(matrix){
	return(apply(matrix, 2, estlambda))
}

qqpercentiles<-function(matrix){
	e = -log10(ppoints(nrow(matrix)))
	mat.sort<-apply(matrix, 2, sort)
	ciValues<--log10(t(apply(mat.sort, 1, quantile,c(0.025,0.5, 0.975))))
	plot(e, ciValues[,2], ylim = c(0, max(ciValues)), type = "n", xlab = "expected", ylab = "observed")
	polygon(c(e, rev(e)), c(ciValues[,1], rev(ciValues[,3])), col = "gray", border = "gray")
	points(e, ciValues[,2], pch = 16, cex = 0.8)
	abline(a = 0, b = 1)
	return(cbind(e, ciValues))
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


## count number of significant associations per simulation
nDMPs<-countNSigMatrix(allSimComb, pThres)
dmpSumStats<-data.frame(matrix(unlist(strsplit(names(nDMPs), "_")), ncol = 2, byrow = TRUE), nDMPs)
colnames(dmpSumStats)[c(1:2)]<-c("Method", "DMPtype")

dmpSumStats$Method[dmpSumStats$DMPtype %in% c("Double-", "NeuN+", "Sox10+")]<-"ctLR"
dmpSumStats$Method<- gsub("MLM", "MER", dmpSumStats$Method)
dmpSumStats$Method<- gsub("LM", "allLR", dmpSumStats$Method)
dmpSumStats$Method<-factor(dmpSumStats$Method, levels = c("ctLR", "allLR", "MER", "CRR"))
dmpSumStats$DMPtype<-factor(dmpSumStats$DMPtype, levels = c("Double-", "NeuN+", "Sox10+", "ME", "Int"))
#boxplot(nDMPs ~ as.factor(names(nDMPs)), ylab = "Number of DMPs", xlab = "Analytical model")
pos <- position_dodge(0.9)

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
  xlab("Regression method")  +
  facet_wrap(vars(DMPtype), nrow = 1, 
  labeller = labeller(DMPtype = c("ME" = "Main effect", "Int" = "Interaction", "Double-" = "DoubleNeg", "NeuN+" = "NeuNPos", "Sox10+" = "Sox10Pos")), scales = "free_x") +
  scale_y_continuous(trans='log10')

ggsave(file.path(dataDir, "Plots", "ViolinPlotnDMPsNullEWAS.pdf"), plot = fig1a, height = 4, width = 12) 


## Simulation Scenario 2: introduce true effects, prespecified to be either cell-type specific or common to all cell-types

outFiles<-list.files(pattern="nSigSimulateTrueEffects_Chunk")
outFiles<-outFiles[grep("pThres", outFiles)]
sumSimComb<-lapply(pThres.opts, function(x) NULL)
for(each in outFiles){
	meanDiff<-unlist(strsplit(strsplit(each, "_MeanDiff")[[1]][2], "\\_pThres_1e-05.rdata"))
	load(each)
	sumSim<-lapply(sumSim, function(mat) {
        cbind(as.numeric(meanDiff), mat)
        colnames(mat)[1]<-"MeanDiff"
        return(mat)
    })
	
	sumSimComb<-mapply(rbind, sumSimComb, sumSim, SIMPLIFY = FALSE)
}

## across simulations calculate mean true positive, false positive and 95% CI per model

sumSimComb<-as.data.frame(sumSimComb)
sumSimComb$propCS<-sumSimComb$nCTspecific/sumSimComb$nProbes


## summarise lambda across proportions of cell-specific associations per model
lambda.sum.propCS<-NULL

for(numDMPs in unique(sumSimComb$nProbes)){
	lambda.sum.propCS[[as.character(numDMPs)]]<-do.call(data.frame, 
	    aggregate(sumSimComb[which(sumSimComb$nProbes == numDMPs),
		c("LM_ME_lambda", "LM_Int_lambda",
		"MLM_ME_lambda", "MLM_Int_lambda",
		"CRR_ME_lambda", "CRR_Int_lambda",
		"LM_Double-_lambda", "LM_NeuN+_lambda", "LM_Sox10+_lambda")], 
	    by = list(sumSimComb$propCS[which(sumSimComb$nProbes == numDMPs)]), 
		FUN = quantile, c(0.025,0.5,0.975)))
}

fig1d<-list()
panelNum<-1
for(numDMPs in unique(sumSimComb$nProbes)){

	tmp<-melt(lambda.sum.propCS[[as.character(numDMPs)]][,c(1,seq(3,19,3))], id = "Group.1")
    tmp<- tmp[grep("_ME_", as.character(tmp$variable)),]
	tmp$variable <- unlist(lapply(strsplit(as.character(tmp$variable), "_"), head, n = 1))
	fig1d[[panelNum]]<-ggplot(tmp, aes(x = Group.1, y = value, colour = variable)) + geom_point(size = 1.5) + geom_line(size = 1.5) + scale_colour_manual(values = methodCols, labels = methodLabels) +ylab("lambda") + xlab("Proportion Cell Specific") + title(paste(numDMPs, "DMPs"))
	panelNum<-panelNum+1
	
	tmp<-melt(lambda.sum.propCS[[as.character(numDMPs)]][,c(1,seq(3,19,3))], id = "Group.1")
    tmp<- tmp[grep("_Int_", as.character(tmp$variable)),]
	tmp$variable <- unlist(lapply(strsplit(as.character(tmp$variable), "_"), head, n = 1))
	fig1d[[panelNum]]<-ggplot(tmp, aes(x = Group.1, y = value, colour = variable)) + geom_point(size = 1.5) + geom_line(size = 1.5) + scale_colour_manual(values = methodCols, labels = methodLabels) +ylab("lambda") + xlab("Proportion Cell Specific")
	panelNum<-panelNum+1

}

fig1d<-ggarrange(plotlist=fig1d[c(seq(1,5,2), seq(2,6,2))], nrow = 2, ncol = 3, common.legend = TRUE, legend.grob = get_legend(fig1a, position = NULL))



y_lim<-range(unlist(lapply(lambda.sum.propCS, range)))
par(mfrow = c(1,3))

	plot(lambda.sum.propCS[[as.character(numDMPs)]][,1], lambda.sum.propCS[[as.character(numDMPs)]][,2][,2], type = "l", ylim = y_lim, 
		xlab = "Proportion DMPs Cell Specific", ylab = "Lambda", col = "red", lty = 1, main = paste0("nDMPs = ", numDMPs))
	lines(lambda.sum.propCS[[as.character(numDMPs)]][,1], lambda.sum.propCS[[as.character(numDMPs)]][,3][,2], col = "red", lty = 2)
	lines(lambda.sum.propCS[[as.character(numDMPs)]][,1], lambda.sum.propCS[[as.character(numDMPs)]][,4][,2], col = "blue", lty = 1)
	lines(lambda.sum.propCS[[as.character(numDMPs)]][,1], lambda.sum.propCS[[as.character(numDMPs)]][,5][,2], col = "blue", lty = 2)
	lines(lambda.sum.propCS[[as.character(numDMPs)]][,1], lambda.sum.propCS[[as.character(numDMPs)]][,6][,2], col = "forestgreen", lty = 1)
	lines(lambda.sum.propCS[[as.character(numDMPs)]][,1], lambda.sum.propCS[[as.character(numDMPs)]][,7][,2], col = "forestgreen", lty = 2)
	abline(h = 1)
}
legend("topright", c("LM_ME", "LM_Int", "MLM_ME", "MLM_Int", "CRR_ME", "CRR_Int"), col = c("red", "red", "blue", "blue", "forestgreen", "forestgreen"), lty = c(1,2))

## summarise lambda as function of total associations per model

lambda.sum.nSig<-aggregate(sumSimComb[,c("LM_ME_lambda", "LM_Int_lambda", "MLM_ME_lambda", "MLM_Int_lambda", "CRR_ME_lambda", "CRR_Int_lambda")], by = list(sumSimComb$nProbes), FUN = quantile, c(0.025,0.5,0.975))

y_lim<-range(lambda.sum.nSig[,-1][,2])

plot(lambda.sum.nSig[,1], lambda.sum.nSig[,2][,2], type = "l", ylim = y_lim, 
	xlab = "nSignificant DMPs", ylab = "Lambda", col = "red", lty = 1)
lines(lambda.sum.nSig[,1], lambda.sum.nSig[,3][,2], col = "red", lty = 2)
lines(lambda.sum.nSig[,1], lambda.sum.nSig[,4][,2], col = "blue", lty = 1)
lines(lambda.sum.nSig[,1], lambda.sum.nSig[,5][,2], col = "blue", lty = 2)
lines(lambda.sum.nSig[,1], lambda.sum.nSig[,6][,2], col = "forestgreen", lty = 1)
lines(lambda.sum.nSig[,1], lambda.sum.nSig[,7][,2], col = "forestgreen", lty = 2)
abline(h = 1)

## summarise lambda as function of total cell-specific associations per model

lambda.sum.nCSSig<-aggregate(sumSimComb[,c("LM_ME_lambda", "LM_Int_lambda", "MLM_ME_lambda", "MLM_Int_lambda", "CRR_ME_lambda", "CRR_Int_lambda")], by = list(sumSimComb$nCTspecific), FUN = quantile, c(0.025,0.5,0.975))

y_lim<-range(lambda.sum.nCSSig[,-1][,2])

plot(lambda.sum.nCSSig[,1], lambda.sum.nCSSig[,2][,2], type = "l", ylim = y_lim, 
	xlab = "nSignificant Cell-specific DMPs", ylab = "Lambda", col = "red", lty = 1)
lines(lambda.sum.nCSSig[,1], lambda.sum.nCSSig[,3][,2], col = "red", lty = 2)
lines(lambda.sum.nCSSig[,1], lambda.sum.nCSSig[,4][,2], col = "blue", lty = 1)
lines(lambda.sum.nCSSig[,1], lambda.sum.nCSSig[,5][,2], col = "blue", lty = 2)
lines(lambda.sum.nCSSig[,1], lambda.sum.nCSSig[,6][,2], col = "forestgreen", lty = 1)
lines(lambda.sum.nCSSig[,1], lambda.sum.nCSSig[,7][,2], col = "forestgreen", lty = 2)
abline(h = 1)

## summarise lambda as function of total cell-consistent associations per model

lambda.sum.nCConSig<-aggregate(sumSimComb[,c("LM_ME_lambda", "LM_Int_lambda", "MLM_ME_lambda", "MLM_Int_lambda", "CRR_ME_lambda", "CRR_Int_lambda")], by = list(sumSimComb$nProbes - sumSimComb$nCTspecific), FUN = quantile, c(0.025,0.5,0.975))

y_lim<-range(lambda.sum.nCConSig[,-1][,2])

plot(lambda.sum.nCConSig[,1], lambda.sum.nCConSig[,2][,2], type = "l", ylim = y_lim, 
	xlab = "nSignificant Cell-consistent DMPs", ylab = "Lambda", col = "red", lty = 1)
lines(lambda.sum.nCConSig[,1], lambda.sum.nCConSig[,3][,2], col = "red", lty = 2)
lines(lambda.sum.nCConSig[,1], lambda.sum.nCConSig[,4][,2], col = "blue", lty = 1)
lines(lambda.sum.nCConSig[,1], lambda.sum.nCConSig[,5][,2], col = "blue", lty = 2)
lines(lambda.sum.nCConSig[,1], lambda.sum.nCConSig[,6][,2], col = "forestgreen", lty = 1)
lines(lambda.sum.nCConSig[,1], lambda.sum.nCConSig[,7][,2], col = "forestgreen", lty = 2)
abline(h = 1)


## compare performance at detecting true associations
## ME analysis should pick up all cell consistent effects and potentially cell specific effects also
## Int analysis should pick up cell-specific effects

sumSimComb <- replace_na(sumSimComb, list(0))
## calculate true positive rate for all DMPs (regardless of type) and regardless of which term in the model
tpRate<-cbind(sumSimComb$MeanDiff, (sumSimComb$LM_ME_nSigTrueDMPs + sumSimComb$LM_Int_nSigTrueDMPs)/(sumSimComb$nProbes), 
		(sumSimComb$MLM_ME_nSigTrueDMPs + sumSimComb$MLM_Int_nSigTrueDMPs)/(sumSimComb$nProbes),
		(sumSimComb$CRR_ME_nSigTrueDMPs + sumSimComb$CRR_Int_nSigTrueDMPs)/(sumSimComb$nProbes),
		sumSimComb$"LM_Double-_nSigTrueDMPs"/(sumSimComb$nProbes - sumSimComb$nCTspecific + sumSimComb$"Double-_nCSTrueDMPs"), 
		sumSimComb$"LM_NeuN+_nSigTrueDMPs"/(sumSimComb$nProbes - sumSimComb$nCTspecific + sumSimComb$"NeuN+_nCSTrueDMPs"), 
		sumSimComb$"LM_Sox10+_nSigTrueDMPs"/(sumSimComb$nProbes - sumSimComb$nCTspecific + sumSimComb$"Sox10+_nCSTrueDMPs"))
colnames(tpRate)<-c("MeanDiff", "LM_LM", "MLM_MLM", "CRR_CRR", "LM_Double-", "LM_NeuN+", "LM_Sox10+")

## calculate true positive rate for consistent effects
tpRate.CCon<-cbind(sumSimComb$MeanDiff, (sumSimComb$LM_ME_nSigTrueCommon)/(sumSimComb$nProbes-sumSimComb$nCTspecific), 
		(sumSimComb$LM_Int_nSigTrueCommon)/(sumSimComb$nProbes-sumSimComb$nCTspecific), 
		(sumSimComb$MLM_ME_nSigTrueCommon)/(sumSimComb$nProbes-sumSimComb$nCTspecific),
		(sumSimComb$MLM_Int_nSigTrueCommon)/(sumSimComb$nProbes-sumSimComb$nCTspecific), 
		(sumSimComb$CRR_ME_nSigTrueCommon)/(sumSimComb$nProbes-sumSimComb$nCTspecific),
		(sumSimComb$CRR_Int_nSigTrueCommon)/(sumSimComb$nProbes-sumSimComb$nCTspecific),
		(sumSimComb$"LM_Double-_nSigTrueCommon")/(sumSimComb$nProbes-sumSimComb$nCTspecific),		
		(sumSimComb$"LM_NeuN+_nSigTrueCommon")/(sumSimComb$nProbes-sumSimComb$nCTspecific),		
		(sumSimComb$"LM_Sox10+_nSigTrueCommon")/(sumSimComb$nProbes-sumSimComb$nCTspecific)		
		)
colnames(tpRate.CCon)<-c("MeanDiff", "LM_ME", "LM_Int", "MLM_ME", "MLM_Int", "CRR_ME", "CRR_Int", "LM_Double-", "LM_NeuN+", "LM_Sox10+")

## calculate true positive rate for ct specific effects
tpRate.CT<-cbind(sumSimComb$MeanDiff, (sumSimComb$LM_ME_nSigTrueCS)/(sumSimComb$nCTspecific),
		(sumSimComb$LM_Int_nSigTrueCS)/(sumSimComb$nCTspecific), 
		(sumSimComb$MLM_ME_nSigTrueCS)/(sumSimComb$nCTspecific),
		(sumSimComb$MLM_Int_nSigTrueCS)/(sumSimComb$nCTspecific),
		(sumSimComb$CRR_ME_nSigTrueCS)/(sumSimComb$nCTspecific),
		(sumSimComb$CRR_Int_nSigTrueCS)/(sumSimComb$nCTspecific),
		(sumSimComb$"Double-_nSigTrueCS")/sumSimComb$"Double-_nCSTrueDMPs", 
		(sumSimComb$"NeuN+_nSigTrueCS")/sumSimComb$"NeuN+_nCSTrueDMPs",
		(sumSimComb$"Sox10+_nSigTrueCS")/sumSimComb$"Sox10+_nCSTrueDMPs"
		)
colnames(tpRate.CT)<-c("MeanDiff", "LM_ME", "LM_Int", "MLM_ME", "MLM_Int", "CRR_ME", "CRR_Int", "LM_Double-", "LM_NeuN+", "LM_Sox10+")

## calculate false positive rate
fpRate<-cbind(sumSimComb$MeanDiff, (sumSimComb$LM_ME_nSigOther)/(sumSimComb$LM_ME_TotSig),
		(sumSimComb$LM_Int_nSigOther)/(sumSimComb$LM_Int_TotSig), 
		(sumSimComb$MLM_ME_nSigOther)/(sumSimComb$MLM_ME_TotSig),
		(sumSimComb$MLM_Int_nSigOther)/(sumSimComb$MLM_Int_TotSig),
		(sumSimComb$CRR_ME_nSigOther)/(sumSimComb$CRR_ME_TotSig),
		(sumSimComb$CRR_Int_nSigOther)/(sumSimComb$CRR_Int_TotSig),
		(sumSimComb$"LM_Double-_nSigOther")/(sumSimComb$"LM_Double-_TotSig"),		
		(sumSimComb$"LM_NeuN+_nSigOther")/(sumSimComb$"LM_NeuN+_TotSig"),		
		(sumSimComb$"LM_Sox10+_nSigOther")/(sumSimComb$"LM_Sox10+_TotSig"))
colnames(fpRate)<-c("MeanDiff", "LM_ME", "LM_Int", "MLM_ME", "MLM_Int", "CRR_ME", "CRR_Int", "LM_Double-", "LM_NeuN+", "LM_Sox10+")

for(meanDiff in c(0.02,0.05)){
	fig2a<-list()
	fig2a[[1]] <- as.data.frame(tpRate) %>% 
	filter(MeanDiff == meanDiff) %>% 
	select(LM_LM:"LM_Sox10+") %>% 
	replace_na(list("LM_Double-" =0, "LM_NeuN+" = 0, "LM_Sox10+" = 0))%>% 
	melt() %>% separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "LM","MLM", "CRR"))) %>% 
	mutate(Regression = factor(Regression)) %>%
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + 
	  scale_fill_manual(values = methodCols, labels = methodLabels, name = "Regression Method")+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) + theme(text = element_text(size = 20))  + ylab("True Positive Rate")    + xlab("Regression method") +labs(title = "All DMPs significant with either term")
				   
	fig2a[[2]] <- as.data.frame(tpRate.CCon) %>% 
	filter(MeanDiff == meanDiff) %>% 
	select(LM_ME:"LM_Sox10+") %>% 
	replace_na(list("LM_Double-" =0, "LM_NeuN+" = 0, "LM_Sox10+" = 0))%>% 
	melt() %>% 
	separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% mutate(Regression = factor(Regression)) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels, name = "Regression Method")+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + ylab("True Positive Rate")    + xlab("Regression method") +labs(title = "Common DMPs")

	fig2a[[3]] <- as.data.frame(tpRate.CT) %>% filter(MeanDiff == meanDiff) %>% select(LM_ME:"LM_Sox10+") %>% replace_na(list("LM_Double-" =0, "LM_NeuN+" = 0, "LM_Sox10+" = 0))%>% melt() %>% separate(variable, c("Regression", "Term")) %>% mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% mutate(Regression = factor(Regression)) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels)+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + 
				   ylab("True Positive Rate")    + xlab("Regression method") +labs(title = "Cell-specific DMPs")
				   
	fig2a[[4]] <- as.data.frame(fpRate)%>% 
	filter(MeanDiff == meanDiff) %>% 
	select(LM_ME:"LM_Sox10+") %>% 
	replace_na(list("LM_Double-" =0, "LM_NeuN+" = 0, "LM_Sox10+" = 0))%>% 
	melt() %>% separate(variable, c("Regression", "Term")) %>% 
	mutate(Term = factor(Term, levels = c("Double", "NeuN", "Sox10", "ME","Int"))) %>% mutate(Regression = factor(Regression)) %>% 
	ggplot( aes(x=Term, y=value, fill = Regression))  +
	  geom_violin(position = pos, scale = 'width')  + scale_fill_manual(values = methodCols, labels = methodLabels)+
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) +  theme(text = element_text(size = 20))  + 
				   ylab("False Positive Rate")    + xlab("Regression method") +labs(title = "")


	fig2a<-ggarrange(plotlist=fig2a, nrow = 2, ncol = 2, common.legend = TRUE)

	ggsave(file.path(dataDir, "Plots", paste0("ViolinPlotCTEWASTPFPRatesMeanDiff", meanDiff, ".pdf")), 
	plot = fig2a, height = 12, width = 12) 

}


fig3a<-list()
tpRate<-as.data.frame(tpRate)
tpRate$nProbes<-sumSimComb$nProbes
tpRate$propCS<-sumSimComb$propCS

## summarise true positive rate as function of total associations per model
fig3a[[1]]<-as.data.frame(group_by(tpRate, nProbes) %>% summarise(LM = mean(LM_LM), MLM = mean(MLM_MLM), CRR = mean(CRR_CRR))) %>% 
melt(id = "nProbes") %>%
ggplot(aes(x = nProbes, y = value, colour = variable)) + 
geom_line(size = 2) +
xlab("Number of DMPs") + 
ylab("True Positive Rate")


## as function of proportion of cell-specific effects
fig3a[[2]]<-as.data.frame(group_by(tpRate, propCS) %>% summarise(LM = mean(LM), MLM = mean(MLM), CRR = mean(CRR))) %>% 
melt(id = "propCS") %>%
ggplot(aes(x = propCS, y = value, colour = variable)) + 
geom_line(size = 2) +
xlab("Proportion CT specific") + 
ylab("True Positive Rate")
fig3a<-ggarrange(plotlist=fig3a, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave(file.path(dataDir, "Plots", "LineGraphCTEWASTPFPRates.pdf"), plot = fig3a, height = 6, width = 12) 

### NOTE number of false positives not rate.
