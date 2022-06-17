##---------------------------------------------------------------------#
##
## Title: Compute Correlation
##
## Purpose of script: to compute the correlation matrix of methylation across sites with >10x coverage, tissue-specific
##
## Author: Jessica Shields
##
## Date Created: 2022-06-08
##
##---------------------------------------------------------------------#

## set working directory
setwd("~/BrainFANS/sequencing")

## clear the R environment
rm(list=ls()) 

## load arguments
args = commandArgs(trailingOnly=TRUE)
project<-args[1]
tissue<-args[2]
sampleids<-args[3:length(args)]

##load config file
#project='WGBS/rizzardi'
source("~/BrainFANS/sequencing/BSSeq/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(data.table)
library(tidyr)
library(ggplot2)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
covFile<-list.files(methylDir, pattern = '.chr1.cov.bg') 
covFile<-covFile[match(sampleids,gsub(".chr1.cov.bg", "", covFile))]


## import chr1 coverage as one dataframe
covMethyl<-NULL
covMethyl<-fread(paste(methylDir, covFile[1], sep = '/'), sep = '\t') %>%
  subset(rowSums(.[,5:6] > 9) > 0)
covMethyl$V7<-sampleids[1]

for (x in 2:length(covFile)){
  tmp<-fread(paste(methylDir, covFile[x], sep = '/'), sep = '\t')
  tmp<-subset(tmp, rowSums(tmp[,5:6] > 9) > 0)
  tmp$V7<-sampleids[x]
  
  #get list of sites common to existing dataframe and sample.dataframe
  commonSites<-merge(covMethyl, tmp, by.x='V2', by.y='V2', all=FALSE)$V2 
  covMethyl<-rbind(covMethyl, tmp)
  
  #filter by common sites
  covMethyl<-covMethyl[which(covMethyl$V2 %in% commonSites)]
} 

## create matrix of samples as cols and sites as rows
covMethyl<- covMethyl[,c(2, 4,7)] %>% #reduce to necessary columns
  pivot_wider(names_from = V7, values_from = V4) 

#----------------------------------------------------------------------#
# GENERATE CORRELATION MATRIX
#----------------------------------------------------------------------#
# number of samples
  ## n(n-1)/2 - formula to calculate number of non-zero elements in n*n triangular matrix
n<-length(colnames(covMethyl))-1

# create empty matrix
corrMat<-NULL
corrMat<- matrix(data = NA, nrow = n, ncol = n-1)
rownames(corrMat)<-colnames(covMethyl)[2:length(colnames(covMethyl))]
colnames(corrMat)<-colnames(covMethyl)[2:(length(colnames(covMethyl))-1)]

# populate matrix
for (y in 1:(n-1)) {
  for (x in (y+1):n) {
    # fill in correlation matrix with site percentages from the imported data  
    corrMat[x,y] <- cor(covMethyl[(x+1)], covMethyl[(y+1)],  method = "pearson", use = "complete.obs")
  }
}

# write to file
write.table(corrMat, paste0(methylDir, '/QCOutput/',tissue,'.corr.qc'), sep='\t', quote = FALSE)


