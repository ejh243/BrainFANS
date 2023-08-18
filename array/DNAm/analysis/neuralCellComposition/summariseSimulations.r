##---------------------------------------------------------------------#
##
## Title: Summarise simulations of  cell deconvolution models
##
## Purpose of script: Aggregate and visualise the simulation results comparing different reference panels 
##
## Author: Eilis Hannon
##
## Date Created: 19/10/2022
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#




#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
rmse<-function(observed, expected){
	sqrt(mean((observed-expected)^2))
}
args<-commandArgs(trailingOnly = TRUE)
resPath<-args[1] # path to simulation results
outFiles<-list.files(pattern = "CompareDeconvolutionParameters", resPath)




#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape)
library(paletteer)


#----------------------------------------------------------------------#
# LOAD AND PLOT DATA
#----------------------------------------------------------------------#

group.colors <- paletteer_d("ggsci::category10_d3")
names(group.colors)<-c("NeuNNeg_SOX10Neg", "NeuNPos", "NeuNNeg", "NeuNNeg_SOX10Pos", "NeuNNeg_Sox10Neg_IRF8Pos", "NeuNNeg_Sox10Neg_IRF8Neg", "SATB2Neg","SATB2Pos", 
"NeuNPos_SOX6Neg", "NeuNPos_SOX6Pos")
panel.colors <- paletteer_d("ggsci::category10_d3", n = 8)
names(panel.colors)<-paste("Panel", 1:8)
group.labels<-c("NeuNNeg/SOX10Neg", "NeuNPos", "NeuNNeg", "NeuNNeg/SOX10Pos", "NeuNNeg/Sox10Neg/IRF8Pos", "NeuNNeg/Sox10Neg/IRF8Neg", "SATB2Neg","SATB2Pos", 
"NeuNPos/SOX6Neg", "NeuNPos/SOX6Pos")

predOutAll<-list()
for(file in outFiles){

	modelNum<-gsub("\\.rdata", "", unlist(strsplit(file, "_"))[2]) 
	load(file.path(resPath, file))
	ct<-unique(colnames(predOut)[duplicated(colnames(predOut))])
	colnames(predOut)[colnames(predOut) %in% ct]<-c(paste0("Prop", ct), paste0("Pred", ct), paste0("Diff", ct))
	predOutAll[[modelNum]]<-predOut
}

fig0a<-list()
fig0b<-list()

cet_lim<-range(unlist(lapply(predOutAll, "[", "CETYGO")))
rmse_lim<-range(unlist(lapply(predOutAll, "[", "RMSE")))

for(i in 1:length(predOutAll)){

	## compare algorithms
	fig0a[[i]]<-ggplot(predOutAll[[i]], aes(x = Selection, y = CETYGO, fill = Selection)) + 
	geom_violin()+ stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")+ 
		xlab("") +
		ylab("CETYGO") + ylim(cet_lim) +
		theme_bw() +  
		theme(legend.position="none")  + ggtitle(paste("Panel", i)) + theme(text = element_text(size = 20))     
		
	fig0b[[i]]<-ggplot(predOutAll[[i]], aes(x = Selection, y = RMSE, fill = Selection)) + 
	geom_violin()+ 
	stat_summary(fun.data=mean_sdl,geom="pointrange", color="black")+ 
		xlab("") +
		ylab("RMSE") + ylim(rmse_lim) +
		theme_bw() +  
		theme(legend.position="none") + ggtitle(paste("Panel", i)) + theme(text = element_text(size = 20))     
}

fig0a[[i]]<- fig0a[[i]] +  
		theme()

fig0b[[i]]<- fig0b[[i]] 

ggarrange(plotlist = fig0a,
			  ncol = 4, nrow = 2, widths = c(4,4,4,4), common.legend = TRUE)
ggsave(file.path(resPath, "plots", paste0("ViolinPlotCETYGOAcrossModels.pdf")), width = 16, height = 8)

## extract sum stats into a table

sumStats<-matrix(data = NA, nrow = length(predOutAll), ncol = 12)
colnames(sumStats)<-c(outer(c("Q25%", "Med", "Q75%"), c("ANOVA_CETYGO", "IDOL_CETYGO", "ANOVA_RMSE", "IDOL_RMSE"), paste))
for(i in 1:length(predOutAll)){
	sumStats[i,]<-c(c(t(aggregate(CETYGO ~ Selection, predOutAll[[i]], quantile, c(0.25,0.5,0.75))$CETYGO)),c(t(aggregate(RMSE ~ Selection, predOutAll[[i]], quantile, c(0.25,0.5,0.75))$RMSE)))
}
write.csv(sumStats, file.path(resPath, paste0("TableSumStatsCETYGORMSEAcrossModels.csv")))

ggarrange(plotlist = fig0b,
			  ncol = 4, nrow = 2, widths = c(4,4,4,4), common.legend = TRUE)
ggsave(file.path(resPath, "plots", paste0("ViolinPlotRMSEAcrossModels.pdf")), width = 16, height = 8)

fig1a<-list()
fig1b<-list()
for(i in 1:length(predOutAll)){
  
	sumCETYGO <- summarySE(predOutAll[[i]], measurevar="CETYGO", groupvars=c("nProbes", "Selection"))
	fig1a[[i]]<-ggplot(sumCETYGO, aes(x=nProbes, y=CETYGO, colour=Selection)) + 
		geom_line() +
		geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
		xlab("Number of sites") +
		ylab("CETYGO") + ylim(0, 0.05) +
		theme_bw() +  
		theme(legend.position="none")+       
		geom_ribbon(aes(ymin=CETYGO-ci, ymax=CETYGO+ci, colour=Selection), linetype=2, alpha=0.1) + ggtitle(paste("Panel", i))

	## compare algorithm parameters
	sumRMSE <- summarySE(predOutAll[[i]], measurevar="RMSE", groupvars=c("nProbes", "Selection"))
	fig1b[[i]]<-ggplot(sumRMSE, aes(x=nProbes, y=RMSE, colour=Selection)) + 
		geom_line() +
		geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
		xlab("Number of sites") +
		ylab("RMSE") + ylim(0,0.15) +
		theme_bw() +  
		theme(legend.position="none")  + 
		geom_ribbon(aes(ymin=RMSE-ci, ymax=RMSE+ci, colour=Selection), linetype=2, alpha=0.1)  + ggtitle(paste("Panel", i))

}
ggarrange(plotlist = fig1a,
			  ncol = 4, nrow = 2, common.legend = TRUE)	
ggsave(file.path(resPath, "plots", "LineGraphCETYGOAgainstnProbesAcrossModels.pdf"), width = 26, height = 15, units = "cm")
ggarrange(plotlist = fig1b,
			  ncol = 4, nrow = 2, common.legend = TRUE)	
ggsave(file.path(resPath, "plots", "LineGraphRMSEAgainstnProbesAcrossModels.pdf"), width = 26, height = 15, units = "cm")

fig3a<-list()
fig3b<-list()
for(i in 1:length(predOutAll)){
  
	predOutByCT<-NULL
	for(each in colnames(predOutAll[[i]])[grep("Prop", colnames(predOutAll[[i]]))]){
		tmp<-predOutAll[[i]][,c("Selection", "nProbes", each, sub("Prop", "Pred", each))]
		tmp$CellType<-gsub("Prop", "", each)
		colnames(tmp)[3:4]<-c("Actual", "Predicted")
		predOutByCT<-rbind(predOutByCT, tmp)
	}

	nCT<-length(grep("Prop", colnames(predOutAll[[i]])))

	summ <- predOutByCT %>% 
	  group_by(Selection,CellType) %>% 
	  summarise(Rsq = cor(Actual, Predicted),
				RMSE = rmse(Actual, Predicted)) %>% 
	  mutate_if(is.numeric, round, digits=2)

	# Here we create our annotations data frame.
	df.annotations <- data.frame()
	# Rsq
	df.annotations <- rbind(df.annotations,
							cbind(as.character(summ$Selection),as.character(summ$CellType),
								  paste("Rsq", summ$Rsq,
										sep = " = ")))
	# RMSE
	df.annotations <- rbind(df.annotations,
							cbind(as.character(summ$Selection),as.character(summ$CellType),
								  paste("RMSE", summ$RMSE,
										sep = " = ")))

	colnames(df.annotations) <- c("Selection", "CellType", "label")

	vertical_adjustment = ifelse(grepl("Rsq",df.annotations$label),1.5,3)


	ggplot(predOutByCT, aes(x = factor(Actual), y = Predicted)) + geom_boxplot() + xlab("Actual") +
	  facet_wrap(~ Selection * CellType, ncol = nCT) + geom_text(data=df.annotations,aes(x=-Inf,y=+Inf,label=label),
				  hjust = -0.1, vjust = vertical_adjustment, size=3.5)  + theme( strip.text = element_text(size = 10))      
				  
	ggsave(file.path(resPath, "plots", paste0("ScatterGraphActualPredicted", i, ".pdf")), width = 6*nCT, height = 15, units = "cm")


	predOutByCT<-NULL
	for(each in colnames(predOutAll[[i]])[grep("Diff", colnames(predOutAll[[i]]))]){
		tmp<-predOutAll[[i]][,c("Selection", "nProbes", each)]
		tmp$CellType<-gsub("Diff", "", each)
		colnames(tmp)[3]<-c("Difference")
		predOutByCT<-rbind(predOutByCT, tmp)
	}


	sumCT <- summarySE(predOutByCT, measurevar="Difference", groupvars=c("nProbes", "Selection", "CellType"))
	y_lim<-c(min(sumCT$Difference-sumCT$ci), max(sumCT$Difference+sumCT$ci))
	fig3a[[i]]<-ggplot(subset(sumCT, Selection == "ANOVA"), aes(x=nProbes, y=Difference, colour=CellType)) + 
			geom_line() +
			geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
			xlab("Number of sites") +
			ylab("Predicted - Actual") +
			ylim(y_lim) +
			theme_bw() +  
		theme(legend.position="none")+       
			geom_ribbon(aes(ymin=Difference-ci, ymax=Difference+ci, colour=CellType), linetype=2, alpha=0.1)  + labs(tag = paste("Panel", i))
	fig3b[[i]]<-ggplot(subset(sumCT, Selection == "IDOL"), aes(x=nProbes, y=Difference, colour=CellType)) + 
			geom_line() +
			geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
			xlab("Number of sites") +
			ylab("Predicted - Actual") +
			ylim(y_lim) +
			theme_bw() +       
			geom_ribbon(aes(ymin=Difference-ci, ymax=Difference+ci, colour=CellType), linetype=2, alpha=0.1)  + labs(tag = paste("Panel", i))
}

ggarrange(plotlist= fig3a,
			  ncol = 4, nrow = 2, widths = c(4,4,4,4), common.legend = TRUE)		
ggsave(file.path(resPath, "plots", "LineGraphCTErrorAgainstnProbesANOVA.pdf"), width = 16, height = 8, units = "cm")
ggarrange(plotlist= fig3b,
			  ncol = 4, nrow = 2, widths = c(4,4,4,4), common.legend = TRUE)		
ggsave(file.path(resPath, "plots", "LineGraphCTErrorAgainstnProbesIDOL.pdf"), width = 16, height = 8, units = "cm")

fig4a<-list()
fig4b<-list()
predOutByCTAll<-NULL
for(i in 1:length(predOutAll)){

	## compare error across cell types		  
	predOutByCT<-melt(predOutAll[[i]], id.vars = c("Selection", "nProbes"),
					measure.vars = colnames(predOutAll[[i]])[grep("Diff", colnames(predOutAll[[i]]))], variable = "CellType", value = "Difference")
	predOutByCT$CellType<-gsub("Diff", "", predOutByCT$CellType)
	
	predOutByCTAll<-rbind(predOutByCTAll, cbind(paste("Panel", i), predOutByCT))

    predOutByCTsub<-subset(predOutByCT, Selection == "ANOVA")
	fig4a[[i]]<-ggplot(predOutByCTsub, aes(x = CellType, y = value, fill = CellType)) + geom_violin() + coord_flip()+ geom_vline(xintercept = 0) + xlab("") 	+ 
	ylim(-0.6, 0.6)  + scale_fill_manual(values=group.colors, labels = group.labels) +  
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.direction="horizontal")+
			ylab("Predicted - Actual") + xlab(paste("Panel", i)) + 
	stat_summary(fun.data=mean_sdl,geom="pointrange", color="black") +
	geom_hline(aes(yintercept = 0))
		
	predOutByCTsub<-subset(predOutByCT, Selection == "IDOL")
	fig4b[[i]]<-ggplot(predOutByCTsub, aes(x = CellType, y = value, fill = CellType)) + geom_violin() + coord_flip()+ geom_vline(xintercept = 0) + xlab("") 	+ 
	ylim(-0.6, 0.6)  + scale_fill_manual(values=group.colors, labels = group.labels) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.direction="horizontal")+
			ylab("Predicted - Actual") + xlab(paste("Panel", i)) + 
	stat_summary(fun.data=mean_sdl,geom="pointrange", color="black") +
	geom_hline(aes(yintercept = 0))
}

ggarrange(plotlist= fig4a,
			  ncol = 4, nrow = 2, common.legend = TRUE, legend.grob = get_legend(fig4a[[1]], position = NULL))	
ggsave(file.path(resPath, "plots", "ViolinplotDifferenceANOVA.pdf"), width = 20, height = 12, units = "cm")

ggarrange(plotlist= fig4b,
			  ncol = 4, nrow = 2, common.legend = TRUE, legend.grob = get_legend(fig4b[[1]], position = NULL))	
ggsave(file.path(resPath, "plots", "ViolinplotDifferenceIDOL.pdf"), width = 20, height = 12, units = "cm")


## plot by cell type
colnames(predOutByCTAll)[1]<-"Panel"
fig5a<-list()
fig5b<-list()
for(ct in names(table(predOutByCTAll$CellType))){

	predOutByCTsub<-subset(predOutByCTAll, Selection == "ANOVA" & CellType == ct)
	fig5a[[ct]]<-ggplot(predOutByCTsub, aes(x = Panel, y = value, fill = Panel)) + geom_violin(width=0.5) + coord_flip()+ geom_vline(xintercept = 0) + xlab("") 	+ 
		ylim(-0.6, 0.6)  +
	  theme(axis.title.y=element_blank(),
			axis.text.y=element_blank(),
			axis.ticks.y=element_blank())+
				ylab("Predicted - Actual") + ggtitle(group.labels[match(ct, names(group.colors))]) + scale_fill_manual(values=panel.colors)  + 
	stat_summary(fun.data=mean_sdl,geom="pointrange", color="black", size = 0.5) +
	geom_hline(aes(yintercept = 0)) 

	predOutByCTsub<-subset(predOutByCTAll, Selection == "IDOL" & CellType == ct)
	fig5b[[ct]]<-ggplot(predOutByCTsub, aes(x = Panel, y = value, fill = Panel)) + geom_violin(width=0.5) + coord_flip()+ geom_vline(xintercept = 0) + xlab("") 	+ 
		ylim(-0.6, 0.6) +
	  theme(axis.title.y=element_blank(),
			axis.text.y=element_blank(),
			axis.ticks.y=element_blank())+
				ylab("Predicted - Actual") + ggtitle(group.labels[match(ct, names(group.colors))]) + scale_fill_manual(values=panel.colors)  + 
	stat_summary(fun.data=mean_sdl,geom="pointrange", color="black", size = 0.5) +
	geom_hline(aes(yintercept = 0))
}

ggarrange(plotlist= fig5a,
			  ncol = 5, nrow = 2, common.legend = TRUE)	
ggsave(file.path(resPath, "plots", "ViolinPlotDifferenceByCellTypeANOVA.pdf"), width = 25, height = 18, units = "cm")

ggarrange(plotlist= fig5b,
			  ncol = 5, nrow = 2, common.legend = TRUE)	
ggsave(file.path(resPath, "plots", "ViolinPlotDifferenceByCellTypeIDOL.pdf"), width = 25, height = 18, units = "cm")


write.csv(aggregate(value ~ CellType * Panel * Selection, predOutByCTAll, summary), file.path(resPath, paste0("TableSumStatsErrorByCelltypePanelSelection.csv")))

write.csv(aggregate(value ~ CellType *Selection, predOutByCTAll, summary), file.path(resPath, paste0("TableSumStatsErrorByCelltypeSelection.csv")))

write.csv(aggregate(value ~ CellType *Selection, subset(predOutByCTAll, !Panel %in% c("Panel 4","Panel 5")), summary), file.path(resPath, paste0("TableSumStatsErrorByCelltypeSelectionNoConflict.csv")))
