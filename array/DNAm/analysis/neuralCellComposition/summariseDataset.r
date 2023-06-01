##---------------------------------------------------------------------#
##
## Title: Summarise dataset
##
## Purpose of script: Create summary tables of demographics and perform PCA analysis 
##
## Author: Eilis Hannon
##
## Date Created: 22/2/2023
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# 

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
## set plotting colours
library(paletteer)
ctCols <- paletteer_d("ggsci::category10_d3")
names(ctCols)<-c("NeuNNeg/SOX10Neg", "NeuNPos", "NeuNNeg", "NeuNNeg/SOX10Pos", "NeuNNeg/Sox10Neg/IRF8Pos", "NeuNNeg/Sox10Neg/IRF8Neg", "SATB2Neg","SATB2Pos", 
"NeuNPos/SOX6Neg", "NeuNPos/SOX6Pos")

args<-commandArgs(trailingOnly = TRUE)
normData <- args[1] # path to R object with normalised datae
plotPath <- args[2]


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)

#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#


load(normData)
norm.all<-norm.all[,pheno.all$Basename]
pos <- position_dodge(0.9)


table(pheno.all$Cell.type)
aggregate(Age ~ Cell.type, FUN = summary, data = pheno.all)
aggregate(Age ~ Cell.type, FUN = sd, data = pheno.all)

## subset to individual level data

indData<-unique(pheno.all[,c("Indidivual.ID", "Age", "Sex", "Tissue.Centre")])
summary(indData)

## PC plots our data only
subCells<-c("NeuNPos", "NeuNNeg/SOX10Pos", "NeuNNeg/SOX10Neg", "NeuNNeg/Sox10Neg/IRF8Pos",  "NeuNNeg/Sox10Neg/IRF8Neg", "SATB2Pos", "SATB2Neg")

norm.sub<-norm.all[,pheno.all$Cell.type %in% subCells]
pheno.sub<-pheno.all[pheno.all$Cell.type %in% subCells,]


pca<-prcomp(t(norm.sub))
pca$sdev[1:10]^2 / sum(pca$sdev^2)

pcDat<-data.frame(pca$x[,1:10], pheno.sub$Cell.type)
colnames(pcDat)[11]<-"Cell.type"

fig1<-list()
fig1[[1]]<-ggplot(pcDat, aes(x=PC1, y = PC2, color = Cell.type)) + geom_point(size = 4) + scale_color_manual(values = ctCols[sort(subCells)],
                           guide = guide_legend(override.aes = list(shape = 15, size = 10) )) + theme(text=element_text(size=21), legend.direction="horizontal") + labs(tag = letters[1]) + labs(color=NULL)


for(i in 1:5){
fig1[[i+1]]<-ggplot(pcDat, aes_string(x="Cell.type", y=paste0("PC",i), fill = "Cell.type"))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white")  + theme(legend.position = "none", text=element_text(size=21), axis.text.x = element_blank()) +
					   scale_fill_manual(values = ctCols[sort(subCells)]) + labs(tag = letters[i+1]) + xlab("")
}


ggarrange(plotlist=fig1, nrow = 2, ncol = 3, common.legend = TRUE, legend.grob = get_legend(fig1[[1]], position = NULL))
ggsave(filename = file.path(plotPath, paste0("PCsByCellFractionExeterOnly.pdf")),  units = "in", width = 20, height = 12)


## PC plots all fractions

pcaAll<-prcomp(t(norm.all))
pcaAll$sdev[1:10]^2 / sum(pcaAll$sdev^2)

pcDatAll<-data.frame(pcaAll$x[,1:10], pheno.all$Cell.type)
colnames(pcDatAll)[11]<-"Cell.type"


fig2<-list()
fig2[[1]]<-ggplot(pcDatAll, aes(x=PC1, y = PC2, color = Cell.type)) + geom_point(size = 4) + scale_color_manual(values = ctCols[sort(subCells)],
                           guide = guide_legend(override.aes = list(shape = 15, size = 10) )) + theme(text=element_text(size=21), legend.direction="horizontal") + labs(tag = letters[1]) + labs(color=NULL)


for(i in 1:5){
fig2[[i+1]]<-ggplot(pcDatAll, aes_string(x="Cell.type", y=paste0("PC",i), fill = "Cell.type"))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white")  + theme(legend.position = "none", text=element_text(size=21), axis.text.x = element_blank()) +
					   scale_fill_manual(values = ctCols[sort(subCells)]) + labs(tag = letters[i+1]) + xlab("")
}


ggarrange(plotlist=fig2, nrow = 2, ncol = 3, common.legend = TRUE, legend.grob = get_legend(fig1[[1]], position = NULL))
ggsave(filename = file.path(plotPath, paste0("PCsByCellFractionCombined.pdf")),  units = "in", width = 18, height = 10)
