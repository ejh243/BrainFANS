##---------------------------------------------------------------------#
##
## Title: Profile cell composition in bulk brain samples
##
## Purpose of script: To calculate cellular composition of bulk brain DNAm
## profiles and characterise against biological factors
##
## Author: Eilis Hannon
##
## Date Created: 13/12/2022
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# takes bulk brain beta matrix and calculates cellular composition
# data provided as R object with betas matrix and pheno data.frame
# requires on execution path to robject with dnam data
# requires on execution path to trained model parameters
# requires on execution path to folder to output plots and tables

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
## set plotting colours
library(paletteer)
ctCols <- paletteer_d("ggsci::category10_d3")
names(ctCols)<-c("DoubleNeg", "NeuNPos", "NEUNNeg", "Sox10Pos", "IRF8Pos", "TripleNeg", "SATB2Neg","SATB2Pos", 
"SOX6Neg", "SOX6Pos")

args<-commandArgs(trailingOnly = TRUE)
bulkPath <- args[1]
modelPath <- args[2]
plotPath <- args[3]

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(CETYGO)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(corrplot)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

load(modelPath)
load(bulkPath)
pos <- position_dodge(0.9)

#----------------------------------------------------------------------#
# SUMMARISE DATASET
#----------------------------------------------------------------------#


if("BrainRegion" %in% colnames(pheno)){
	sumStat<-cbind(table(pheno$BrainRegion),
		aggregate(Age ~ BrainRegion, pheno, mean)$Age,
		aggregate(Age ~ BrainRegion, pheno, sd)$Age,
		aggregate(Age ~ BrainRegion, pheno, min)$Age,
		aggregate(Age ~ BrainRegion, pheno, max)$Age,
		table(pheno$BrainRegion, pheno$Sex))

	write.csv(sumStat, file = file.path(plotPath, "SummariseDataset.csv"))
}

if("AgeBin" %in% colnames(pheno)){
	sumStat<-cbind(table(droplevels(pheno$AgeBin)), 
	aggregate(Age ~ AgeBin, pheno,min)$Age, 
	aggregate(Age ~ AgeBin, pheno,max)$Age, 
	table(droplevels(pheno$AgeBin),pheno$Sex))
		write.csv(sumStat, file = file.path(plotPath, "SummariseDatasetByAgeBin.csv"))
}

#----------------------------------------------------------------------#
# CALCULATE CELL COMPOSITION
#----------------------------------------------------------------------#

predCCANOVA<-list()
predCCIDOL<-list()
for(i in 1:length(brainCoefANOVA)){
	sites<-intersect(rownames(betas), rownames(brainCoefANOVA[[i]]))
	predCCANOVA[[i]]<-projectCellTypeWithError(betas[sites,], brainCoefANOVA[[i]][sites,])
	if(!is.null(brainCoefIDOL[[i]])){
		sites<-intersect(rownames(betas), rownames(brainCoefIDOL[[i]]))
		predCCIDOL[[i]]<-projectCellTypeWithError(betas[sites,], brainCoefIDOL[[i]][sites,])
	}
}
save(predCCANOVA, predCCIDOL, file = paste0(strsplit(bulkPath, "\\.")[[1]][1], "CellCompEstimates.Rdata"))

## reformat into single data.frame for plotting
sumOut<-NULL
for(i in 1:length(predCCANOVA)){
	sumOut<-rbind(sumOut, data.frame("Method" = "ANOVA", "Model" = i, "BrainRegion" = pheno$BrainRegion, "Age" = pheno$Age, "Sex" = pheno$Sex, "Batch" = pheno$Sentrix_ID, "CETYGO" = predCCANOVA[[i]][,"CETYGO"], "SumProp" = rowSums(predCCANOVA[[i]][,1:(ncol(predCCANOVA[[i]])-2)])))
	if(!is.null(predCCIDOL[[i]])){
		sumOut<-rbind(sumOut, data.frame("Method" = "IDOL", "Model" = i, "BrainRegion" = pheno$BrainRegion, "Age" = pheno$Age, "Sex" = pheno$Sex, "Batch" = pheno$Sentrix_ID, "CETYGO" = predCCIDOL[[i]][,"CETYGO"], "SumProp" = rowSums(predCCIDOL[[i]][,1:(ncol(predCCIDOL[[i]])-2)])))
	}
}

sumOut$Model<-as.factor(sumOut$Model)

#----------------------------------------------------------------------#
# TEST CETYGO AGAINST BRAIN REGION, SEX & AGE
#----------------------------------------------------------------------#

pheno$Age2<-pheno$Age^2

## reset PFC as baseline region
if(length(unique(pheno$BrainRegion)) > 1){  
	lmOut<-NULL
	for(i in 1:length(predCCANOVA)){
		model<-lm(predCCANOVA[[i]][,"CETYGO"] ~ Age + Age2 + Sex + relevel(as.factor(pheno$BrainRegion), "BA9"), data = pheno)
		lmOut<-rbind(lmOut, c(i, c(t(summary(model)$coefficients[-1,c(1,4)]))))
	}
	colnames(lmOut)<-c("Panel","AgeEst", "Age2P", "Age2Est", "Age2P", "SexEst", "SexP", c(outer(c("Est_", "P_"), levels(relevel(as.factor(pheno$BrainRegion), "BA9"))[-1], paste0)))
	write.csv(lmOut, file = paste0(plotPath, "RegressionAgainstCETYGOANOVA.csv"))
	lmOut<-NULL
	for(i in 1:length(predCCIDOL)){
		if(!is.null(predCCIDOL[[i]])){
		model<-lm(predCCIDOL[[i]][,"CETYGO"] ~ Age + Age2 + Sex + relevel(as.factor(pheno$BrainRegion), "BA9"), data = pheno)
		lmOut<-rbind(lmOut, c(i, c(t(summary(model)$coefficients[-1,c(1,4)]))))
		}
	}
	colnames(lmOut)<-c("Panel","AgeEst", "Age2P", "Age2Est", "Age2P", "SexEst", "SexP", c(outer(c("Est_", "P_"), levels(relevel(as.factor(pheno$BrainRegion), "BA9"))[-1], paste0)))
	write.csv(lmOut, file = paste0(plotPath, "RegressionAgainstCETYGOIDOL.csv"))
	
	## exclude CER
	lmOut<-NULL
	for(i in 1:length(predCCANOVA)){
		model<-lm(predCCANOVA[[i]][,"CETYGO"] ~ Age + Age2 + Sex + relevel(as.factor(pheno$BrainRegion), "BA9"), data = pheno, subset = pheno$BrainRegion != "CEREB")
		lmOut<-rbind(lmOut, c(i, c(t(summary(model)$coefficients[-1,c(1,4)]))))
	}
	colnames(lmOut)<-c("Panel","AgeEst", "Age2P", "Age2Est", "Age2P", "SexEst", "SexP", c(outer(c("Est_", "P_"), levels(droplevels(relevel(as.factor(pheno$BrainRegion), "BA9"), exclude = "CEREB"))[-1], paste0)))
	write.csv(lmOut, file = paste0(plotPath, "RegressionAgainstCETYGOANOVAExcludeCER.csv"))
	lmOut<-NULL
	for(i in 1:length(predCCIDOL)){
		if(!is.null(predCCIDOL[[i]])){
		model<-lm(predCCIDOL[[i]][,"CETYGO"] ~ Age + Age2 + Sex + relevel(as.factor(pheno$BrainRegion), "BA9"), data = pheno, subset = pheno$BrainRegion != "CEREB")
		lmOut<-rbind(lmOut, c(i, c(t(summary(model)$coefficients[-1,c(1,4)]))))
		}
	}
	colnames(lmOut)<-c("Panel","AgeEst", "Age2P", "Age2Est", "Age2P", "SexEst", "SexP", c(outer(c("Est_", "P_"), levels(droplevels(relevel(as.factor(pheno$BrainRegion), "BA9"), exclude = "CEREB"))[-1], paste0)))
	write.csv(lmOut, file = paste0(plotPath, "RegressionAgainstCETYGOIDOLExcludeCER.csv"))
	
} else {
	lmOut<-NULL
	for(i in 1:length(predCCANOVA)){
		model<-lm(predCCANOVA[[i]][,"CETYGO"] ~ Age + Age2 + Sex + Sentrix_ID, data = pheno)
		null<-lm(predCCANOVA[[i]][,"CETYGO"] ~ Age + Age2 + Sex, data = pheno)
		lmOut<-rbind(lmOut, c(i, c(t(summary(model)$coefficients[c("Age", "Age2", "SexM"),c(1,4)]), anova(model, null)[2,6])))
	}
	colnames(lmOut)<-c("Panel","AgeEst", "AgeP", "Age2Est", "Age2P", "SexEst", "SexP", "BatchP")
	write.csv(lmOut, file = paste0(plotPath, "RegressionAgainstCETYGOANOVA.csv"))
	
	lmOut<-NULL
	for(i in 1:length(predCCIDOL)){
		if(!is.null(predCCIDOL[[i]])){
			model<-lm(predCCIDOL[[i]][,"CETYGO"] ~ Age + Age2 + Sex + Sentrix_ID, data = pheno)
			null<-lm(predCCIDOL[[i]][,"CETYGO"] ~ Age + Age2 + Sex, data = pheno)
			lmOut<-rbind(lmOut, c(i, c(t(summary(model)$coefficients[c("Age", "Age2", "SexM"),c(1,4)]), anova(model, null)[2,6])))
		}
	}
	colnames(lmOut)<-c("Panel","AgeEst", "AgeP", "Age2Est", "Age2P", "SexEst", "SexP", "BatchP")
	write.csv(lmOut, file = paste0(plotPath, "RegressionAgainstCETYGOIDOL.csv"))
}

#----------------------------------------------------------------------#
# PLOT CETYGO AGAINST BRAIN REGION
#----------------------------------------------------------------------#

y_lim<-range(sumOut$CETYGO)
fig2a <- ggplot(na.omit(subset(sumOut, Method == "ANOVA")), aes(x=Model, y=CETYGO, fill = BrainRegion))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
  ylim(y_lim)
  
 fig2b <- ggplot(na.omit(subset(sumOut, Method == "IDOL")), aes(x=Model, y=CETYGO, fill = BrainRegion))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
  ylim(y_lim)
ggarrange(fig2a, fig2b, nrow = 2, ncol = 1)
ggsave(filename = file.path(plotPath, "ViolinPlotCETYGOAcrossPanels.pdf"),  units = "in", width = 12, height = 8)

for(br in unique(pheno$BrainRegion)){
	fig3a<-ggplot(na.omit(subset(sumOut, BrainRegion == br & Method == "ANOVA")), aes(x=Model, y=CETYGO))  +
		geom_violin(position = pos, scale = 'width')  +
		stat_summary(fun = "mean", 
		   geom = "point", 
		   position = pos) + 
		geom_vline(xintercept = c(1:length(predCCANOVA))+ 0.5, linetype="dotted") +
		ylim(y_lim) + theme(legend.position = "none")
	fig3b<-ggplot(na.omit(subset(sumOut, BrainRegion == br & Method == "IDOL")), aes(x=Model, y=CETYGO))  +
		geom_violin(position = pos, scale = 'width')  +
		stat_summary(fun = "mean", 
		   geom = "point", 
		   position = pos) + 
		geom_vline(xintercept = c(1:length(predCCANOVA))+ 0.5, linetype="dotted") +
		ylim(y_lim) + theme(legend.position = "none")
	ggarrange(plotlist=list(fig3a, fig3b), nrow = 2)
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotCETYGOAcrossPanels",br,"Samples.pdf")),  units = "in", width = 12, height = 8)
}

#----------------------------------------------------------------------#
# PLOT CETYGO AGAINST AGE
#----------------------------------------------------------------------#

fig5a<-list()
fig5b<-list()
for(i in 1:length(predCCANOVA)){
	tmp_long<-as.data.frame(predCCANOVA[[i]])
	tmp_long$Age<-pheno$Age

	fig5a[[i]] <- ggplot(tmp_long, aes(x=Age, y=CETYGO))  +
	  geom_point() +
     geom_smooth(method = "lm", alpha = .15) + theme(legend.position = "none")

	if(!is.null(predCCIDOL[[i]])){
	tmp_long<-as.data.frame(predCCIDOL[[i]])
	tmp_long$Age<-pheno$Age

	fig5b[[i]] <- ggplot(tmp_long, aes(x=Age, y=CETYGO))  +
	  geom_point() +
     geom_smooth(method = "lm", alpha = .15) + theme(legend.position = "none")
	}else{
			fig5b[[i]]<-ggplot(tmp_long, aes(x=Age, y=CETYGO)) + geom_blank() + theme(legend.position = "none")
	}
}

ggarrange(plotlist=fig5a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossPanelsANOVA.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig5b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossPanelsIDOL.pdf"),  units = "in", width = 12, height = 8)


if("AgeBin" %in% colnames(pheno)){
	fig2a <- ggplot(na.omit(subset(sumOut, Method == "ANOVA")), aes(x=Model, y=CETYGO, fill = AgeBin))  +
	  geom_violin(position = pos, scale = 'width')  +
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
	  ylim(y_lim)
	  
	 fig2b <- ggplot(na.omit(subset(sumOut, Method == "IDOL")), aes(x=Model, y=CETYGO, fill = AgeBin))  +
	  geom_violin(position = pos, scale = 'width')  +
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
	  ylim(y_lim)
	ggarrange(fig2a, fig2b, nrow = 2, ncol = 1)
	ggsave(filename = file.path(plotPath, "ViolinPlotCETYGOAcrossPanelsByAgeBin.pdf"),  units = "in", width = 12, height = 8)
}

#----------------------------------------------------------------------#
# PLOT CETYGO AGAINST SEX
#----------------------------------------------------------------------#

y_lim<-range(sumOut$CETYGO)
fig2a <- ggplot(na.omit(subset(sumOut, Method == "ANOVA")), aes(x=Model, y=CETYGO, fill = Sex))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
  ylim(y_lim)
  
 fig2b <- ggplot(na.omit(subset(sumOut, Method == "IDOL")), aes(x=Model, y=CETYGO, fill = Sex))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
  ylim(y_lim)
ggarrange(fig2a, fig2b, nrow = 2, ncol = 1)
ggsave(filename = file.path(plotPath, "ViolinPlotCETYGOAcrossPanelsBySex.pdf"),  units = "in", width = 12, height = 8)

#----------------------------------------------------------------------#
# PLOT CETYGO AGAINST BATCH
#----------------------------------------------------------------------#

y_lim<-range(sumOut$CETYGO)
fig2a <- ggplot(na.omit(subset(sumOut, Method == "ANOVA")), aes(x=Model, y=CETYGO, fill = Batch))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
  ylim(y_lim)
  
 fig2b <- ggplot(na.omit(subset(sumOut, Method == "IDOL")), aes(x=Model, y=CETYGO, fill = Batch))  +
  geom_violin(position = pos, scale = 'width')  +
  stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + geom_vline(xintercept = c(1:length(predCCANOVA))+0.5, linetype="dotted") +
  ylim(y_lim)
ggarrange(fig2a, fig2b, nrow = 2, ncol = 1)
ggsave(filename = file.path(plotPath, "ViolinPlotCETYGOAcrossPanelsByBatch.pdf"),  units = "in", width = 12, height = 8)


#----------------------------------------------------------------------#
# PLOT DISTRIBUTION OF PREDICTED CELLULAR COMPOSITION: PANEL 8
#----------------------------------------------------------------------#

# calc neuronal & glial proportions

relCC<-predCCIDOL[[8]][,c("IRF8Pos", "SOX6Neg", "SOX6Pos", "Sox10Pos","TripleNeg")]
relCC<-cbind(relCC[,"SOX6Neg"]+relCC[,"SOX6Pos"], relCC[,"IRF8Pos"]+relCC[,"Sox10Pos"]+relCC[,"TripleNeg"], relCC)
colnames(relCC)[1:2]<-c("neuronal", "glial")

tmp<-data.frame(predCCIDOL[[8]])
tmp_long <- gather(tmp, "PredCT", "Proportion")
tmp_long <- subset(tmp_long, PredCT %in% names(ctCols))
tmp_long$BrainRegion<-pheno$BrainRegion
tmp_long<-na.omit(tmp_long)
ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols[unique(tmp_long$PredCT)], drop = TRUE)

ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnSamplesPanel8.pdf")),  units = "in", width = 18, height = 8)


ctMeans<-cbind(colMeans(relCC), apply(relCC, 2, sd))

write.csv(ctMeans, file = file.path(plotPath, "CellCompositionPanel8IDOLSummaryStats.csv"))
#----------------------------------------------------------------------#
# CORPLOT OF ALL CELL TYPE PREDICTIONS 
#----------------------------------------------------------------------#


matPredANOVA<-NULL
matPredIDOL<-NULL
matError<-NULL
for(i in 1:length(predCCANOVA)){
	nCT<-ncol(predCCANOVA[[i]])-2
	colnames(predCCANOVA[[i]])[1:nCT] <- paste0("M", i, "_", colnames(predCCANOVA[[i]])[1:nCT])
	matPredANOVA<-cbind(matPredANOVA, predCCANOVA[[i]][,1:nCT])
	matError<-cbind(matError, predCCANOVA[[i]][,"CETYGO"])
	colnames(matError)[ncol(matError)]<-paste0("ANOVA_M", i)
	if(!is.null(predCCIDOL[[i]])){
		colnames(predCCIDOL[[i]])[1:nCT] <- paste0("IDOL_M", i, "_", colnames(predCCIDOL[[i]])[1:nCT])
		matPredIDOL<-cbind(matPredIDOL, predCCIDOL[[i]][,1:nCT])
		matError<-cbind(matError, predCCIDOL[[i]][,"CETYGO"])
		colnames(matError)[ncol(matError)]<-paste0("IDOL_M", i)
	}
}

corMatPred<-cor(matPredANOVA)
pdf(file.path(plotPath, "CorMatPredCompositionANOVA.pdf"), width = 15, height =12)
corrplot.mixed(corMatPred, order = 'hclust',tl.pos = "lt")
dev.off()


corMatPred<-cor(matPredIDOL)
pdf(file.path(plotPath, "CorMatPredCompositionIDOL.pdf"), width = 15, height =12)
corrplot.mixed(corMatPred, order = 'hclust',tl.pos = "lt")
dev.off()


corMatError<-cor(matError)
pdf(file.path(plotPath, "CorMatCETYGO.pdf"), width = 15, height =12)
corrplot.mixed(corMatError, order = 'hclust',tl.pos = "lt")
dev.off()




#----------------------------------------------------------------------#
# EXTRACT BEST PREDICTION FOR EACH CELL TYPE
#----------------------------------------------------------------------#

predCCBest<-cbind(predCCANOVA[[1]][,c("M1_DoubleNeg", "M1_NeuNPos")], predCCIDOL[[5]][,c("IDOL_M5_IRF8Pos","IDOL_M5_TripleNeg")], predCCANOVA[[3]][,c("M3_SATB2Neg", "M3_SATB2Pos")], predCCANOVA[[6]][,c("M6_NEUNNeg", "M6_SOX6Pos", "M6_SOX6Neg")], predCCIDOL[[4]][,"IDOL_M4_Sox10Pos"])
colnames(predCCBest)[10]<-"Sox10Pos"
colnames(predCCBest)<-gsub("M._", "", colnames(predCCBest))
colnames(predCCBest)<-gsub("IDOL_", "", colnames(predCCBest))


ctMeans<-cbind(colMeans(predCCBest), apply(predCCBest, 2, sd))

write.csv(ctMeans, file = file.path(plotPath, "CellCompositionBestSummaryStats.csv"))
#----------------------------------------------------------------------#
# TEST AGAINST AGE,SEX, BRAIN REGION
#----------------------------------------------------------------------#

lmOut<-NULL
if(length(unique(pheno$BrainRegion)) > 1){  
	for(j in 1:ncol(predCCBest)){
		model<-lm(predCCBest[,j] ~ Age + Age2 + Sex + relevel(as.factor(pheno$BrainRegion), "BA9"), data = pheno)
		lmOut<-rbind(lmOut, c(colnames(predCCBest)[j], c(t(summary(model)$coefficients[-1,c(1,4)]))))
	}
	colnames(lmOut)<-c("Panel","AgeEst", "Age2P", "Age2Est", "Age2P", "SexEst", "SexP", c(outer(c("Est_", "P_"), levels(relevel(as.factor(pheno$BrainRegion), "BA9"))[-1], paste0)))
	write.csv(lmOut, file = file.path(plotPath, "RegressionofAgeSexBrainRegionAgainstCellCompositionBest.csv"))
} else {
	for(j in 1:ncol(predCCBest)){
		model<-lm(predCCBest[,j] ~ Age + Age2 + Sex + Sentrix_ID, data = pheno)
		null<-lm(predCCBest[,j] ~ Age + Age2 + Sex, data = pheno)
		lmOut<-rbind(lmOut, c(colnames(predCCBest)[j], c(t(summary(model)$coefficients[c("Age", "Age2", "SexM"),c(1,2,4)])), anova(model,null)[2,6]))
	}

	colnames(lmOut)<-c("CT", "AgeEst", "AgeSE", "AgeP", "Age2Est", "Age2SE", "Age2P", "SexEst", "SexSE", "SexP", "BatchP")
	write.csv(lmOut, file = file.path(plotPath, "RegressionofAgeSexAgainstCellCompositionBest.csv"))
}

#----------------------------------------------------------------------#
# PLOT DISTRIBUTION OF PREDICTED CELLULAR COMPOSITION
#----------------------------------------------------------------------#
 
tmp<-data.frame(predCCBest)
tmp_long<- gather(tmp, "PredCT", "Proportion")
tmp_long$BrainRegion<-pheno$BrainRegion
tmp_long<-na.omit(tmp_long)
ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols)

ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnSamplesBest.pdf")),  units = "in", width = 18, height = 8)

ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols) +
			facet_wrap(vars(BrainRegion), nrow = 2)
			
ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnSamplesBestByRegion.pdf")),  units = "in", width = 18, height = 8)

tmp_long$Sex<-pheno$Sex
ggplot(tmp_long, aes(x=Sex, y=Proportion, fill = Sex))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + 
			facet_wrap(vars(PredCT), nrow = 2)
			
ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnSamplesBestBySex.pdf")),  units = "in", width = 18, height = 8)

for(br in unique(pheno$BrainRegion)){
	tmp<-data.frame(predCCBest[which(pheno$BrainRegion == br), ])
	tmp_long<- gather(tmp, "PredCT", "Proportion")
	tmp_long$Age<-pheno$Age[which(pheno$BrainRegion == br)]
	tmp_long$Sex<-pheno$Sex[which(pheno$BrainRegion == br)]
	ggplot(tmp_long, aes(x=Age, y=Proportion, color = PredCT))  +
			geom_point() +
			geom_smooth(method = "lm", alpha = .15, aes(fill = PredCT)) + theme(legend.position = "none") + scale_color_manual(values = ctCols) +
	  facet_wrap(vars(PredCT), nrow = 2)

	ggsave(filename = file.path(plotPath, paste0("ScatterplotPredPropnAgainstAgeSamplesBest",br,".pdf")),  units = "in", width = 18, height = 8)
	
	ggplot(tmp_long, aes(x=Sex, y=Proportion, fill = Sex))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + 
			facet_wrap(vars(PredCT), nrow = 2)
			
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnSamplesBestBySex",br,".pdf")),  units = "in", width = 18, height = 8)
}

