##---------------------------------------------------------------------#
##
## Title: Compute Correlation Part 2
##
## Purpose of script: to compute the correlation matrix of methylation across sites with >10x coverage, tissue-specific
##
## Author: Jessica Shields
##
## Date Created: 2022-06-08
##
##---------------------------------------------------------------------#

# ## set working directory
# setwd("~/BrainFANS/sequencing")

## clear the R environment
rm(list=ls()) 

## load arguments
args = commandArgs(trailingOnly=TRUE)
project<-args[1]
tissue<-args[2]
sampleids<-args[3:length(args)]

print(c(project, tissue, sampleids))

##load config file (if not submitted to isca, set example args)
if (length(args)==0){
  project<-'WGBS/epiGaba'
  tissue<- 'dorsolateral_prefrontal_cortex'
  sampleids<- c("SRR5343780",
                "SRR5343781",
                "SRR5343788",
                "SRR5343795",
                "SRR5343802",
                "SRR5343803",
                "SRR5343810",
                "SRR5343817",
                "SRR5343818",
                "SRR5343830",
                "SRR5343838",
                "SRR5343840",
                "SRR5343845",
                "SRR5343848")
  sampleids<- c('GABA1_BS',
                'GABA2_BS',
                'GLU1_BS',
                'GLU2_BS',
                'OLIG1_BS',
                'OLIG2_BS')
  tissue<- 'anterior_cingulate_cortex_BA24'
  sampleids<- c("GABA1_BS",
                "GABA2_BS",
                "GLU1_BS",
                "GLU2_BS",
                "OLIG1_BS",
                "OLIG2_BS")
}

source("BSSeq/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(data.table)
library(tidyr)
library(ggplot2)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

covFile<-list.files(qcDir, pattern = '.chr1.cov.bg') 

print(sampleids)
print(gsub(".chr1.cov.bg", "", covFile))
print(match(sampleids,gsub(".chr1.cov.bg", "", covFile)))
covFile<-covFile[match(sampleids,gsub(".chr1.cov.bg", "", covFile))]


## import chr1 coverage as one dataframe
covMethyl<-NULL
covMethyl<-fread(paste(qcDir, covFile[1], sep = '/'), sep = '\t') %>%
  subset(rowSums(.[,5:6] > 9) > 0)
covMethyl$V7<-sampleids[1]

for (x in 2:length(covFile)){
  tmp<-fread(paste(qcDir, covFile[x], sep = '/'), sep = '\t')
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
covMethyl<-covMethyl[,order(colnames(covMethyl[1:ncol(covMethyl)]))]

# create empty matrix
corrMat<-NULL
corrMat<- matrix(data = NA, nrow = n, ncol = n-1)
rownames(corrMat)<-colnames(covMethyl)[1:length(colnames(covMethyl))-1]
colnames(corrMat)<-colnames(covMethyl)[1:(length(colnames(covMethyl))-2)]


# populate matrix
for (y in 1:(n-1)) {
  for (x in (y+1):n) {
    print(colnames(covMethyl[x]))
    print(colnames(covMethyl[y]))
    # fill in correlation matrix with site percentages from the imported data  
    corrMat[colnames(covMethyl[x]),colnames(covMethyl[y])] <- cor(covMethyl[x], covMethyl[y],  method = "pearson", use = "complete.obs")
  }
}

# write to file
write.table(corrMat, paste0(methylDir, '/QCOutput/',tissue,'.corr.qc'), sep='\t', quote = FALSE)


