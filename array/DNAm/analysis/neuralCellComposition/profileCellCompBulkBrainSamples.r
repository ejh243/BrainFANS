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

# remove fetal samples from brain region analysis
pheno$BrainRegion[which(pheno$Age < 0)]<-NA


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
	sumOut<-rbind(sumOut, data.frame("Method" = "ANOVA", "Model" = i, "BrainRegion" = pheno$BrainRegion, "Age" = pheno$Age, "CETYGO" = predCCANOVA[[i]][,"CETYGO"], "SumProp" = rowSums(predCCANOVA[[i]][,1:(ncol(predCCANOVA[[i]])-2)])))
	if(!is.null(predCCIDOL[[i]])){
		sumOut<-rbind(sumOut, data.frame("Method" = "IDOL", "Model" = i, "BrainRegion" = pheno$BrainRegion, "Age" = pheno$Age, "CETYGO" = predCCIDOL[[i]][,"CETYGO"], "SumProp" = rowSums(predCCIDOL[[i]][,1:(ncol(predCCIDOL[[i]])-2)])))
	}
}

sumOut$Model<-as.factor(sumOut$Model)

#----------------------------------------------------------------------#
# CREATE PLOTS
#----------------------------------------------------------------------#
pos <- position_dodge(0.9)

## For each brain region plot distribution of predicted cellular composition
for(br in unique(pheno$BrainRegion)){
	fig1a<-list()
	for(i in 1:length(predCCANOVA)){
		cellTypes <- colnames(predCCANOVA[[i]])[1:(ncol(predCCANOVA[[i]])-2)]
		tmp<-data.frame(predCCANOVA[[i]][which(pheno$BrainRegion == br), cellTypes])
		tmp_long<- gather(tmp, "PredCT", "Proportion", cellTypes)
		fig1a[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols[cellTypes]) + theme(legend.position = "none")
	}
	ggarrange(plotlist=fig1a, nrow = 2, ncol = 4)
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnAcrossModels",br,"SamplesANOVA.pdf")),  units = "in", width = 18, height = 8)

	fig1b<-list()
	for(i in 1:length(predCCIDOL)){
		if(!is.null(predCCIDOL[[i]])){
		cellTypes <- colnames(predCCIDOL[[i]])[1:(ncol(predCCIDOL[[i]])-2)]
		tmp<-data.frame(predCCIDOL[[i]][which(pheno$BrainRegion == br), cellTypes])
		tmp_long<- gather(tmp, "PredCT", "Proportion", cellTypes)
		fig1b[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols[cellTypes]) + theme(legend.position = "none")
		} else{
			fig1b[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT)) + geom_blank()
		}
	}
	ggarrange(plotlist=fig1b, nrow = 2, ncol = 4)
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnAcrossModels",br,"SamplesIDOL.pdf")),  units = "in", width = 18, height = 8)
}

## plot for each panel all brain regions 

for(i in 1:length(predCCANOVA)){
	fig1a<-list()
	cellTypes <- colnames(predCCANOVA[[i]])[1:(ncol(predCCANOVA[[i]])-2)]
	for(br in unique(pheno$BrainRegion)){
		tmp<-data.frame(predCCANOVA[[i]][which(pheno$BrainRegion == br), cellTypes])
		tmp_long<- gather(tmp, "PredCT", "Proportion", cellTypes)
		fig1a[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols[cellTypes]) + theme(legend.position = "none")
	}
	ggarrange(plotlist=fig1a, nrow = 2, ncol = 4)
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnAcrossRegionsPanel",i,"SamplesANOVA.pdf")),  units = "in", width = 18, height = 8)

	fig1b<-list()
	for(i in 1:length(predCCIDOL)){
		if(!is.null(predCCIDOL[[i]])){
		cellTypes <- colnames(predCCIDOL[[i]])[1:(ncol(predCCIDOL[[i]])-2)]
		tmp<-data.frame(predCCIDOL[[i]][which(pheno$BrainRegion == br), cellTypes])
		tmp_long<- gather(tmp, "PredCT", "Proportion", cellTypes)
		fig1b[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols[cellTypes]) + theme(legend.position = "none")
		} else{
			fig1b[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT)) + geom_blank()
		}
	}
	ggarrange(plotlist=fig1b, nrow = 2, ncol = 4)
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotPredPropnAcrossRegionsPanel",i,"SamplesIDOL.pdf")),  units = "in", width = 18, height = 8)
}


## plot CETYGO grouped by brain region

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
ggsave(filename = file.path(plotPath, "ViolinPlotCETYGOAcrossModels.pdf"),  units = "in", width = 12, height = 8)

## plot CETYGO separately for each brain region

for(br in unique(pheno$BrainRegion)){
	fig3a<-list()
	for(i in 1:length(predCCANOVA)){
		pheno$BrainRegion == br
		fig3a[[i]]<-ggplot(na.omit(subset(sumOut, BrainRegion == br)), aes(x=Model, y=CETYGO, fill = Method))  +
			geom_violin(position = pos, scale = 'width')  +
			stat_summary(fun = "mean", 
               geom = "point", 
               position = pos) + 
			geom_vline(xintercept = c(1:length(predCCANOVA))+ 0.5, linetype="dotted") +
			ylim(y_lim) + theme(legend.position = "none")
	}
	ggarrange(plotlist=fig3a, nrow = 2, ncol = 4)
	ggsave(filename = file.path(plotPath, paste0("ViolinPlotCETYGOAcrossModels",br,"Samples.pdf")),  units = "in", width = 12, height = 8)

}

## plot cell composition against age
fig4a<-list()
fig4b<-list()
for(i in 1:length(predCCANOVA)){
	cellTypes <- colnames(predCCANOVA[[i]])[1:(ncol(predCCANOVA[[i]])-2)]
	tmp_long<-gather(as.data.frame(predCCANOVA[[i]]), "PredCT", "Proportion", cellTypes)
	tmp_long$Age<-pheno$Age
	
	fig4a[[i]] <- ggplot(tmp_long, aes(x=Age, y=Proportion, color = PredCT))  +
	  geom_point() +
     geom_smooth(method = "lm", alpha = .15, aes(fill = PredCT)) + theme(legend.position = "none") + scale_color_manual(values = ctCols[cellTypes])


	if(!is.null(predCCIDOL[[i]])){
		tmp_long<-gather(as.data.frame(predCCIDOL[[i]]), "PredCT", "Proportion", cellTypes)
		tmp_long$Age<-pheno$Age	
		fig4b[[i]] <- ggplot(tmp_long, aes(x=Age, y=Proportion, color = PredCT))  +
		geom_point() +
		geom_smooth(method = "lm", alpha = .15, aes(fill = PredCT)) + theme(legend.position = "none") + scale_color_manual(values = ctCols[cellTypes])
	}else{
			fig4b[[i]]<-ggplot(tmp_long, aes(x=PredCT, y=Proportion, fill = PredCT)) + geom_blank() + theme(legend.position = "none")
	}
	  
 }
ggarrange(plotlist=fig4a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsANOVA.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig4b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsIDOL.pdf"),  units = "in", width = 12, height = 8)

## zoom in on childhood

for(i in 1:length(predCCANOVA)){
	fig4a[[i]] <- fig4a[[i]] + xlim(-1,18)
	if(!is.null(predCCIDOL[[i]])){
		fig4b[[i]] <- fig4b[[i]] + xlim(-1,18)
	}
}
 ggarrange(plotlist=fig4a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsANOVAUnder18.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig4b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsIDOLUnder18.pdf"),  units = "in", width = 12, height = 8)

## fetal and early postnatal

for(i in 1:length(predCCANOVA)){
	fig4a[[i]] <- fig4a[[i]] + xlim(-0.75,2)
	if(!is.null(predCCIDOL[[i]])){
		fig4b[[i]] <- fig4b[[i]] + xlim(-0.75,2)
	}
}
 ggarrange(plotlist=fig4a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsANOVAUnder2.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig4b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsIDOLUnder2.pdf"),  units = "in", width = 12, height = 8)

## omit childhood

for(i in 1:length(predCCANOVA)){
	fig4a[[i]] <- fig4a[[i]] + xlim(18, 110)
	if(!is.null(predCCIDOL[[i]])){
		fig4b[[i]] <- fig4b[[i]] + xlim(18, 110)
	}
}
 ggarrange(plotlist=fig4a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsANOVAAdult.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig4b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotPredPropnagainstAgeAcrossModelsIDOLAdult.pdf"),  units = "in", width = 12, height = 8)

## plot error against age

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
## zoom in on childhood
 ggarrange(plotlist=fig5a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsANOVA.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig5b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsIDOL.pdf"),  units = "in", width = 12, height = 8)
for(i in 1:length(predCCANOVA)){
	fig5a[[i]] <- fig5a[[i]] + xlim(-1,18)
	fig5b[[i]] <- fig5b[[i]] + xlim(-1,18)
}
 ggarrange(plotlist=fig5a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsANOVAUnder18.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig5b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsIDOLUnder18.pdf"),  units = "in", width = 12, height = 8)

## omit childhood

for(i in 1:length(predCCANOVA)){
	fig5a[[i]] <- fig5a[[i]] + xlim(18, 110)
	fig5b[[i]] <- fig5b[[i]] + xlim(18, 110)
}
 ggarrange(plotlist=fig5a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsANOVAAdult.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig5b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsIDOLAdult.pdf"),  units = "in", width = 12, height = 8)


## fetal and early postnatal

for(i in 1:length(predCCANOVA)){
	fig5a[[i]] <- fig5a[[i]] + xlim(-0.75,2)
	fig5b[[i]] <- fig5b[[i]] + xlim(-0.75,2)
}

ggarrange(plotlist=fig4a, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsANOVAUnder2.pdf"),  units = "in", width = 12, height = 8)

ggarrange(plotlist=fig4b, nrow = 2, ncol = 4)
ggsave(filename = file.path(plotPath, "ScatterPlotCETYGOagainstAgeAcrossModelsIDOLUnder2.pdf"),  units = "in", width = 12, height = 8)


## corplot of all cell type predictions 
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
pdf(file.path(plotPath, "CorMatPredCompositionANOVA.pdf"), width = 12, height =12)
corrplot.mixed(corMatPred, order = 'hclust',tl.pos = "lt")
dev.off()


corMatPred<-cor(matPredIDOL)
pdf(file.path(plotPath, "CorMatPredCompositionIDOL.pdf"), width = 12, height =12)
corrplot.mixed(corMatPred, order = 'hclust',tl.pos = "lt")
dev.off()


corMatError<-cor(matError)
pdf(file.path(plotPath, "CorMatCETYGO.pdf"), width = 12, height =12)
corrplot.mixed(corMatError, order = 'hclust',tl.pos = "lt")
dev.off()

## limit to adults

corMatPred<-cor(matPredANOVA[which(pheno$Age > 18),])
pdf(file.path(plotPath, "CorMatPredCompositionANOVAAdults.pdf"), width = 12, height =12)
corrplot.mixed(corMatPred, order = 'hclust',tl.pos = "lt")
dev.off()


corMatPred<-cor(matPredIDOL[which(pheno$Age > 18),])
pdf(file.path(plotPath, "CorMatPredCompositionIDOLAdults.pdf"), width = 12, height =12)
corrplot.mixed(corMatPred, order = 'hclust',tl.pos = "lt")
dev.off()


corMatError<-cor(matError[which(pheno$Age > 18),])
pdf(file.path(plotPath, "CorMatCETYGOAdults.pdf"), width = 12, height =12)
corrplot.mixed(corMatError, order = 'hclust',tl.pos = "lt")
dev.off()