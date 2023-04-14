##---------------------------------------------------------------------#
##
## Title: Summarise dataset
##
## Purpose of script: 
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
names(ctCols)<-c("DoubleNeg", "NeuNPos", "NEUNNeg", "Sox10Pos", "IRF8Pos", "TripleNeg", "SATB2Neg","SATB2Pos", 
"SOX6Neg", "SOX6Pos")

args<-commandArgs(trailingOnly = TRUE)
normData <- args[1]
sampleFile <- args[2]
plotPath <- args[3]


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)

#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#


load(normData)
pos <- position_dodge(0.9)


pheno.all<-pheno.all[which(pheno.all$Cell.type != "Total"),]
norm.all<-norm.all[,pheno.all$Basename]

pheno.all$Cell.type<-gsub("\\+", "Pos", pheno.all$Cell.type) ## need to remove the "+" and "-"
pheno.all$Cell.type<-gsub("\\-", "Neg", pheno.all$Cell.type) ## need to remove the "+" and "-"

sampleSheet<-read.csv(sampleFile)
if(!"Basename" %in% colnames(sampleSheet)){
	sampleSheet$Basename<-paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
}

table(pheno.all$Cell.type)
aggregate(Age ~ Cell.type, FUN = summary, data = sampleSheet)
aggregate(Age ~ Cell.type, FUN = sd, data = sampleSheet)

## subset to individual level data

indData<-unique(sampleSheet[,c("Indidivual.ID", "Age", "Sex", "Tissue.Centre")])

## PC plots our data only
subCells<-c("NeuNPos", "Sox10Pos", "DoubleNeg", "IRF8Pos",  "TripleNeg", "SATB2Pos", "SATB2Neg")

norm.sub<-norm.all[,pheno.all$Cell.type %in% subCells]
pheno.sub<-pheno.all[pheno.all$Cell.type %in% subCells,]


pca<-prcomp(t(norm.sub))
pca$sdev[1:10]^2 / sum(pca$sdev^2)

pcDat<-data.frame(pca$x[,1:10], pheno.sub$Cell.type)
colnames(pcDat)[11]<-"Cell.type"

fig1<-list()
fig1[[1]]<-ggplot(pcDat, aes(x=PC1, y = PC2, color = Cell.type)) + geom_point() + scale_color_manual(values = ctCols[sort(subCells)]) + theme(text=element_text(size=21)) + labs(tag = letters[1])


for(i in 1:5){
fig1[[i+1]]<-ggplot(pcDat, aes_string(x="Cell.type", y=paste0("PC",i), fill = "Cell.type"))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white")  + theme(legend.position = "none", text=element_text(size=21), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
					   scale_fill_manual(values = ctCols[sort(subCells)]) + labs(tag = letters[i+1]) + xlab("")
}


ggarrange(plotlist=fig1, nrow = 2, ncol = 3)
ggsave(filename = file.path(plotPath, paste0("PCsByCellFractionExeterOnly.pdf")),  units = "in", width = 20, height = 12)


## PC plots all fractions

pcaAll<-prcomp(t(norm.all))
pcaAll$sdev[1:10]^2 / sum(pcaAll$sdev^2)

pcDatAll<-data.frame(pcaAll$x[,1:10], pheno.all$Cell.type)
colnames(pcDatAll)[11]<-"Cell.type"

fig2<-list()
fig2[[1]]<-ggplot(pcDatAll, aes(x=PC1, y = PC2, color = Cell.type)) + geom_point() + scale_color_manual(values = ctCols) + theme(text=element_text(size=21))+ labs(tag = letters[1])


for(i in 1:5){
fig2[[i+1]]<-ggplot(pcDatAll, aes_string(x="Cell.type", y=paste0("PC",i), fill = "Cell.type"))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white")  + theme(legend.position = "none", text=element_text(size=21), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
					   scale_fill_manual(values = ctCols) + labs(tag = letters[i+1]) + xlab("")
}


ggarrange(plotlist=fig2, nrow = 2, ncol = 3)
ggsave(filename = file.path(plotPath, paste0("PCsByCellFractionCombined.pdf")),  units = "in", width = 18, height = 10)
