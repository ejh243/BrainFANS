##---------------------------------------------------------------------#
##
## Title: Plot variance explained statistics 
##
## Purpose of script: To summarise variance explained by cellular 
## composition.
##
## Author: Eilis Hannon
##
## Date Created: 03/02/2023
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
## set plotting colours
library(paletteer)
ctCols <- paletteer_d("ggsci::category10_d3")
names(ctCols)<-c("DoubleNeg", "NeuNPos", "NEUNNeg", "Sox10Pos", "IRF8Pos", "TripleNeg", "SATB2Neg","SATB2Pos", 
"SOX6Neg", "SOX6Pos")


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)

resPath <- args[1]
outName<- args[2]


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(variancePartition)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

load(file.path(resPath, "varExplained", paste0(outName,  "VarExplainedCellComp.Rdata")))

# calculate total variance explained by all cell types
for(i in 1:length(varPartANOVA)){
    varPartANOVA[[i]]<-as.data.frame(varPartANOVA[[i]])
    varPartANOVA[[i]]$SumCellComp<-rowSums(varPartANOVA[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
	
	if(!is.null(varPartIDOL[[i]])){
		varPartIDOL[[i]]<-as.data.frame(varPartIDOL[[i]])
		varPartIDOL[[i]]$SumCellComp<-rowSums(varPartIDOL[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
	}
	
	pcaPartANOVA[[i]]<-as.data.frame(pcaPartANOVA[[i]])
    pcaPartANOVA[[i]]$SumCellComp<-rowSums(pcaPartANOVA[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
	
	if(!is.null(pcaPartIDOL[[i]])){	
		pcaPartIDOL[[i]]<-as.data.frame(pcaPartIDOL[[i]])
		pcaPartIDOL[[i]]$SumCellComp<-rowSums(pcaPartIDOL[[i]] %>% select(-c("Age", "Sex", "Sentrix_ID", "Residuals", "CETYGO")))
	}
}


#----------------------------------------------------------------------#
# PLOT RESULTS
#----------------------------------------------------------------------#
pos <- position_dodge(0.9)

fig1a<-list()
fig2a<-list()
fig3a<-list()
fig4a<-list()
fig1b<-list()
fig2b<-list()
fig3b<-list()
fig4b<-list()
allSum<-list()
for(i in 1:length(varPartANOVA)){
	# plot distribution grouped by factors
	tmp<-gather(varPartANOVA[[i]] %>% select(c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals")))
	fig1a[[i]]<-ggplot(tmp, aes(x=key, y=value, fill = key))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "Proportion Var. Explained")
	
	# plot distribution of each cell type
	tmp<-gather(varPartANOVA[[i]] %>% select(-c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals")))
	cellTypes<-unique(tmp$key)
	fig2a[[i]]<-ggplot(tmp, aes(x=key, y=value, fill = key))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + scale_fill_manual(values = ctCols[cellTypes]) + theme(legend.position = "none") + labs(x = "", y = "Proportion Var. Explained")
	
	allSum[["ANOVA"]]<-cbind(allSum[["ANOVA"]], varPartANOVA[[i]]$SumCellComp)
	
	fig3a[[i]]<-rownames_to_column(pcaPartANOVA[[i]]) %>% 
	select(c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals", "rowname")) %>% 
	gather( key = "Factor", value = "Proportion", -rowname) %>%	
	ggplot( aes(x = rowname, y = Proportion, fill = Factor)) +
	  geom_col() +
	  scale_fill_brewer(palette = "Set2") +
	  theme_minimal(base_size = 16) +
	  ylab("Percentage") +
	  xlab(NULL) + theme(legend.position = "none")
	  
	fig4a[[i]]<-rownames_to_column(pcaPartANOVA[[i]]) %>% 
	select(-c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals")) %>% 
	gather( key = "CellType", value = "Proportion", -rowname) %>%	
	ggplot( aes(x = rowname, y = Proportion, fill = CellType)) +
	  geom_col() + scale_fill_manual(values = ctCols[cellTypes]) +
	  theme_minimal(base_size = 16) +
	  ylab("Percentage") +
	  xlab(NULL) + theme(legend.position = "none")
	  
	if(!is.null(varPartIDOL[[i]])){
		tmp<-gather(varPartIDOL[[i]] %>% select(c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals")))
		fig1b[[i]]<-ggplot(tmp, aes(x=key, y=value, fill = key))  +
			  geom_violin(position = pos, scale = 'width')  +
			  stat_summary(fun = "mean", 
						   geom = "point", 
						   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "Proportion Var. Explained")
		
		# plot distribution of each cell type
		tmp<-gather(varPartIDOL[[i]] %>% select(-c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals")))
		cellTypes<-unique(tmp$key)
		fig2b[[i]]<-ggplot(tmp, aes(x=key, y=value, fill = key))  +
			  geom_violin(position = pos, scale = 'width')  +
			  stat_summary(fun = "mean", 
						   geom = "point", 
						   position = pos, col = "white") + scale_fill_manual(values = ctCols[cellTypes]) + theme(legend.position = "none") + labs(x = "", y = "Proportion Var. Explained")
		
		allSum[["IDOL"]]<-cbind(allSum[["IDOL"]], varPartIDOL[[i]]$SumCellComp)
		
		fig3b[[i]]<-rownames_to_column(pcaPartIDOL[[i]]) %>% 
		select(c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals", "rowname")) %>% 
		gather( key = "Factor", value = "Proportion", -rowname) %>%	
		ggplot( aes(x = rowname, y = Proportion, fill = Factor)) +
		  geom_col() +
		  scale_fill_brewer(palette = "Set2") +
		  theme_minimal(base_size = 16) +
		  ylab("Percentage") +
		  xlab(NULL) + theme(legend.position = "none")
		  
		fig4b[[i]]<-rownames_to_column(pcaPartIDOL[[i]]) %>% 
		select(-c("SumCellComp", "Age", "Sex", "Sentrix_ID", "CETYGO", "Residuals")) %>% 
		gather( key = "CellType", value = "Proportion", -rowname) %>%	
		ggplot( aes(x = rowname, y = Proportion, fill = CellType)) +
		  geom_col() + scale_fill_manual(values = ctCols[cellTypes]) +
		  theme_minimal(base_size = 16) +
		  ylab("Percentage") +
		  xlab(NULL) + theme(legend.position = "none")
	}
}

ggarrange(plotlist=fig1a, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "ViolinPlotPropExpldAllFactorsANOVA.pdf")),  units = "in", width = 18, height = 8)

ggarrange(plotlist=fig2a, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "ViolinPlotPropExpldCellTypesANOVA.pdf")),  units = "in", width = 18, height = 8)
	
ggarrange(plotlist=fig3a, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "StackedBarplotPropExpldAllFactorsANOVA.pdf")),  units = "in", width = 18, height = 8)
	
ggarrange(plotlist=fig4a, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "StackedBarplotPropExpldCellTypesANOVA.pdf")),  units = "in", width = 18, height = 8)


ggarrange(plotlist=fig1b, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "ViolinPlotPropExpldAllFactorsIDOL.pdf")),  units = "in", width = 18, height = 8)

ggarrange(plotlist=fig2b, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "ViolinPlotPropExpldCellTypesIDOL.pdf")),  units = "in", width = 18, height = 8)
	
ggarrange(plotlist=fig3b, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "StackedBarplotPropExpldAllFactorsIDOL.pdf")),  units = "in", width = 18, height = 8)
	
ggarrange(plotlist=fig4b, nrow = 2, ncol = 4)
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "StackedBarplotPropExpldCellTypesIDOL.pdf")),  units = "in", width = 18, height = 8)

colnames(allSum[["ANOVA"]])<-paste("Panel", 1:8)
gather(as.data.frame(allSum[["ANOVA"]])) %>%
ggplot(aes(x=key, y=value, fill = key))  +
			  geom_violin(position = pos, scale = 'width')  +
			  stat_summary(fun = "mean", 
						   geom = "point", 
						   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "Proportion Var. Explained")
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "ViolinPropExpldAllCellTypesacrossPanelsANOVA.pdf")),  units = "in", width = 18, height = 8)


colnames(allSum[["IDOL"]])<-paste("Panel", c(1:2,4:8))
gather(as.data.frame(allSum[["IDOL"]])) %>%
ggplot(aes(x=key, y=value, fill = key))  +
			  geom_violin(position = pos, scale = 'width')  +
			  stat_summary(fun = "mean", 
						   geom = "point", 
						   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "Proportion Var. Explained")
ggsave(filename = file.path(resPath, "varExplained", "plots", paste0(outName,  "ViolinPropExpldAllCellTypesacrossPanelsIDOL.pdf")),  units = "in", width = 18, height = 8)
