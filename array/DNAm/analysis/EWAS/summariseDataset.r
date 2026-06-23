##---------------------------------------------------------------------#
##
## Title: Summarise dataset
##
## Purpose of script: Create summary tables of demographics for samples included 
## in case control analyses
##
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
normData <- args[1] # path to R object with normalised data


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)

pos <- position_dodge(0.9)

#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#

load(normData)

## remove cell types with less than 20 samples
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]


#----------------------------------------------------------------------#
# CREATE DEMOGRAPHICS TABLE PER SAMPLE TYPE
#----------------------------------------------------------------------#


table(QCmetrics$Cell.type)
table(QCmetrics$Phenotype == "Schizophrenia")
table(QCmetrics$Phenotype == "Schizophrenia", QCmetrics$Cell.type)
table(QCmetrics$Sex == "M")
table(QCmetrics$Sex == "M", QCmetrics$Cell.type)
table(QCmetrics$Sex == "M", QCmetrics$Phenotype)
chisq.test(table(QCmetrics$Sex == "M", QCmetrics$Phenotype))
table(QCmetrics$Sex, QCmetrics$Phenotype, QCmetrics$Cell.type)
lapply(split(QCmetrics, QCmetrics$Cell.type), function(subset){
    chisq.test(table(subset$Phenotype, subset$Sex))$p.value
    }
)

mean(QCmetrics$Age, na.rm = TRUE)
sd(QCmetrics$Age, na.rm = TRUE)

aggregate(Age ~ Cell.type, FUN = mean, data = QCmetrics, na.rm = TRUE)
aggregate(Age ~ Cell.type, FUN = sd, data = QCmetrics, na.rm = TRUE)

aggregate(Age ~ Phenotype, FUN = mean, data = QCmetrics, na.rm = TRUE)
aggregate(Age ~ Phenotype, FUN = sd, data = QCmetrics, na.rm = TRUE)

t.test(Age ~ Phenotype, data = QCmetrics, na.rm = TRUE)

aggregate(Age ~ Cell.type*Phenotype, FUN = mean, data = QCmetrics, na.rm = TRUE)
aggregate(Age ~ Cell.type*Phenotype, FUN = sd, data = QCmetrics, na.rm = TRUE)

lapply(split(QCmetrics, QCmetrics$Cell.type), function(subset){
    t.test(Age ~ Phenotype, data = subset)$p.value
    }
)

#----------------------------------------------------------------------#
# PLOT DISTRIBUTIONS
#----------------------------------------------------------------------#

datLong<-QCmetrics[,c("Cell.type", "Phenotype", "Age")]
mean_data <- datLong %>%
  group_by(Cell.type, Phenotype) %>%
  summarize(mean_Age = mean(Age, na.rm = TRUE))
ggplot(datLong, aes(Phenotype, Age, fill = Phenotype))+
          geom_violin(scale = "width", position= position_dodge(width = 0.9))+ facet_wrap (~ Cell.type)  + 
		   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = mean_data, aes(y = mean_Age, group = Phenotype),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.9)) 
ggsave(file.path(resPath, "Plots", "ViolinPlotAgeByCellTypePhenotype.pdf"),
       width = 7, height = 5, dpi = 150, units = "in")


#----------------------------------------------------------------------#
# CREATE DEMOGRAPHICS TABLE AT INDIVIDUAL LEVEL
#----------------------------------------------------------------------#

## subset to individual level data
indData<-unique(QCmetrics[,c("Indidivual.ID", "Age", "Sex", "Tissue.Centre", "Phenotype")])
summary(indData)
