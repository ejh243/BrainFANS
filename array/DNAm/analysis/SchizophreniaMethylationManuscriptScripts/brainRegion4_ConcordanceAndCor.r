##---------------------------------------------------------------------#
##
## Title: Brain region correlations 
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# correlate neun EWAS results between all PFC and subset PFC for hippocampus and striatum

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#


library(ggplot2)
library(dplyr)
library(ggpubr)

dataDir <- "Schizophrenia/brainRegionData/"
resDir <- file.path(dataDir, "4_analysis/EWAS")
outPath <- file.path(dataDir, "4_analysis/plots")

normData <- file.path(dataDir, "3_normalised/normalised.rdata")



# load results
allRes <- read.csv(file.path(dataDir, "/4_analysis/plots/NeuN+FullResultsAllBrainRegionsAllProbes.csv"))

discovery <- allRes[allRes$P.fullPFC < 1e-04,]


# define function to test concordance
testConcordance <- function(df, col1, col2){
    res <- list()
    df$concordant <- NA
    for(i in 1:nrow(df)){
        if(df[,col1][i] > 0 & df[,col2][i] > 0){
            df$concordant[i] <- "yes"
        } else if(df[,col1][i] < 0 & df[,col2][i] < 0){
            df$concordant[i] <- "yes"
        } else {
            df$concordant[i] <- "no"
        }
    }


    res[[1]] <- table(df$concordant)
    res[[2]] <- prop.table(table(df$concordant))
    res[[3]] <- cor(df[,col1], df[,col2])
    return(res)
}



# hipp concordant

# full PFC
test <- discovery[,c("ProbeID", "Coeff.fullPFC", "Coeff.Hipp")]
test <- test[complete.cases(test), ]
dim(test)

testConcordance(test, "Coeff.fullPFC", "Coeff.Hipp")



# striatum concordant

# 1e-04 full PFC
test <- discovery[,c("ProbeID", "Coeff.fullPFC", "Coeff.Striatum")]
test <- test[complete.cases(test), ]
dim(test)

testConcordance(test, "Coeff.fullPFC", "Coeff.Striatum")



# matched PFC vs full PFC

test <- discovery[,c("ProbeID", "Coeff.matchedPFC", "Coeff.fullPFC")]
test <- test[complete.cases(test), ]
dim(test)

testConcordance(test, "Coeff.matchedPFC", "Coeff.fullPFC")


