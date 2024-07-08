##---------------------------------------------------------------------#
##
## Title: Test cell type composition against disease status
##
## Purpose of script: calculate cellular composition of brain cell types
## and test for differences between cases and controls
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE CETYGO FUNCTION FOR BRAIN TISSUE
#----------------------------------------------------------------------#

calcBrainComposition <- function(betas){
  
  predPropAll<-list() # store output in list

  # run CETYGO using both ANOVA and IDOL methods and all available cell ref panels
    for(method in names(modelBrainCoef)){
        for(j in 1:length(modelBrainCoef[[method]])){
        if(!is.null(modelBrainCoef[[method]][[j]])){
                predPropAll[[method]][[j]]<-projectCellTypeWithError(betas, modelBrainCoef[[method]][[j]])
        }
        } 
    }
    predPropBest<-cbind(predPropAll[["ANOVA"]][[1]][,c("NeuNNeg_SOX10Neg", "NeuNPos")], predPropAll[["IDOL"]][[5]][,c("NeuNNeg_Sox10Neg_IRF8Pos","NeuNNeg_Sox10Neg_IRF8Neg")], predPropAll[["ANOVA"]][[3]][,c("SATB2Neg", "SATB2Pos")], predPropAll[["ANOVA"]][[6]][,c("NeuNNeg", "NeuNPos_SOX6Pos", "NeuNPos_SOX6Neg")], predPropAll[["IDOL"]][[4]][,"NeuNNeg_SOX10Pos"])
    colnames(predPropBest)[10]<-"NeuNNeg_SOX10Pos"

    # add in cetygo scores
    predPropBest<-as.data.frame(predPropBest)
    predPropBest$CETYGO<-predPropAll[["ANOVA"]][[1]][,c("CETYGO")]
    return(predPropBest)
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(CETYGO)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(dplyr)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
cellTypes<-unique(QCmetrics$Cell.type)

predPropBest<-calcBrainComposition(celltypeNormbeta)

save(predPropBest, file = file.path(resPath, "EstimatedNeuralCellComposition.rdata"))

#----------------------------------------------------------------------#
# PLOT DISTRIBUTIONS
#----------------------------------------------------------------------#

predPropBest<-cbind(QCmetrics, predPropBest)

datLong<-gather(predPropBest[,c("Cell.type", "Phenotype", "NeuNNeg_SOX10Neg","NeuNPos","NeuNNeg_Sox10Neg_IRF8Pos","NeuNNeg_Sox10Neg_IRF8Neg","SATB2Neg","SATB2Pos","NeuNNeg","NeuNPos_SOX6Pos","NeuNPos_SOX6Neg","NeuNNeg_SOX10Pos", "CETYGO")], Fraction, Proportion, NeuNNeg_SOX10Neg:CETYGO)

sumStats<-cbind(aggregate(Proportion ~ ., data = datLong, FUN = mean),
aggregate(Proportion ~ ., data = datLong, FUN = sd)[,4])
write.csv(sumStats, file.path(resPath, "EWAS", "Tables", "CellCompositionSummaryStatisticsByCellTypePhenotype.csv"))

testPhenotype<-NULL
for(fraction in c("NeuNNeg_SOX10Neg","NeuNPos","NeuNNeg_Sox10Neg_IRF8Pos","NeuNNeg_Sox10Neg_IRF8Neg","SATB2Neg","SATB2Pos","NeuNNeg","NeuNPos_SOX6Pos","NeuNPos_SOX6Neg","NeuNNeg_SOX10Pos", "CETYGO")){
	for(sampletype in unique(predPropBest$Cell.type)){
		subData<-subset(datLong, Cell.type == sampletype & Fraction == fraction)
		testDat<-t.test(Proportion ~ Phenotype, data = subData)
		testPhenotype<-rbind(testPhenotype, c(sampletype, fraction, testDat$estimate, testDat$p.value))
	}

}
colnames(testPhenotype)<-c("Cell.type", "Fraction", "MeanControls", "MeanSchizophrenia", "Pvalue")
write.csv(testPhenotype, file.path(resPath,  "EWAS", "Tables", "CellCompositionTtestStatisticsByCellTypePhenotype.csv"))

# confirm purity of each cell type
datLong<-gather(predPropBest[,c("Cell.type", "Phenotype", "NeuNPos", "NeuNNeg_SOX10Pos", "NeuNNeg_SOX10Neg")], Fraction, Proportion, NeuNPos:NeuNNeg_SOX10Neg)
mean_data <- datLong %>%
  group_by(Cell.type, Fraction, Phenotype) %>%
  summarize(mean_Proportion = mean(Proportion))
ggplot(datLong, aes(Fraction, Proportion, fill = Phenotype))+
          geom_violin(scale = "width", position= position_dodge(width = 0.9)) + facet_wrap (~ Cell.type) +
		   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = mean_data, aes(y = mean_Proportion, group = Phenotype),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = seq(0, max(datLong$Proportion), by = 0.2))
ggsave(file.path(resPath,  "EWAS", "Plots", "ViolinPlotCellCompositionByCellTypePhenotype.pdf"),
       width = 15, height = 5, dpi = 150, units = "in")


# Plot CETYGO by phenotype

datLong<-predPropBest[,c("Cell.type", "Phenotype", "CETYGO")]
mean_data <- datLong %>%
  group_by(Cell.type, Phenotype) %>%
  summarize(mean_Proportion = mean(CETYGO))
ggplot(datLong, aes(Phenotype, CETYGO, fill = Phenotype))+
          geom_violin(scale = "width", position= position_dodge(width = 0.9))+ facet_wrap (~ Cell.type)  + 
		   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = mean_data, aes(y = mean_Proportion, group = Phenotype),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.9)) 
ggsave(file.path(resPath, "EWAS", "Plots", "ViolinPlotCETYGOByCellTypePhenotype.pdf"),
       width = 7, height = 5, dpi = 150, units = "in")


## Within Double negative look at IRF8 Pos vs Neg
datLong<-gather(predPropBest[,c("Cell.type", "Phenotype", "NeuNNeg_Sox10Neg_IRF8Pos","NeuNNeg_Sox10Neg_IRF8Neg")], Fraction, Proportion, NeuNNeg_Sox10Neg_IRF8Pos:NeuNNeg_Sox10Neg_IRF8Neg) %>% subset(Cell.type == "Double-")
mean_data <- datLong %>%
  group_by(Fraction, Phenotype) %>%
  summarize(mean_Proportion = mean(Proportion))
ggplot(datLong, aes(Fraction, Proportion, fill = Phenotype))+
          geom_violin(scale = "width", position= position_dodge(width = 0.9))+ 
		   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = mean_data, aes(y = mean_Proportion, group = Phenotype),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = seq(0, max(datLong$Proportion), by = 0.2))
ggsave(file.path(resPath,  "EWAS", "Plots", "ViolinPlotCellCompositionByCellTypePhenotypeGlialSubtypes.pdf"),
       width = 5, height = 5, dpi = 150, units = "in")

## Within Neurons look at SOX6 Pos vs Neg
datLong<-gather(predPropBest[,c("Cell.type", "Phenotype", "NeuNPos_SOX6Pos","NeuNPos_SOX6Neg")], Fraction, Proportion, NeuNPos_SOX6Pos:NeuNPos_SOX6Neg) %>% subset(Cell.type == "NeuN+")
mean_data <- datLong %>%
  group_by(Fraction, Phenotype) %>%
  summarize(mean_Proportion = mean(Proportion))
ggplot(datLong, aes(Fraction, Proportion, fill = Phenotype))+
          geom_violin(scale = "width", position= position_dodge(width = 0.9))+ 
		   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = mean_data, aes(y = mean_Proportion, group = Phenotype),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = seq(0, max(datLong$Proportion), by = 0.2))
ggsave(file.path(resPath, "EWAS", "Plots",  "ViolinPlotCellCompositionByCellTypePhenotypeNeuronalSubtypes.pdf"),
       width = 5, height = 5, dpi = 150, units = "in")