##---------------------------------------------------------------------#
##
## Title: Simulate chromHMM data =======================================
##
## Purpose of script: to simulate binary matrix for chromHMM input 
##
## Author: Jessica Shields 
##
## Date Created: 2022-12-13
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 
## set working directory
setwd("/lustre/home/jms260/BrainFANS")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args<-"roadMapPFC"
} 
intproject<-args
simproject<-'simulateRM'
type <- 'PFC' #'N+'

source("integrative/chromHMM/config/config.r")

# create pathway for simulated data
simulateDir<-gsub("[^/]+$", paste0(simproject, "/2_mergeBinarised"), chromDir) # catch everything after the last "/"

#dir.create(simulateDir, recursive=TRUE)

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#
library(magrittr)


readEmissionTable<-function(tablePath){
  t<-read.table(tablePath, sep='\t', header = TRUE, row.names = 1) %>% 
    as.matrix() %>%
    #.[, colOrder]# %>% # reorder to same col order
    {(.>=0.5)*1} # binarise dataframe (threshold 0.5)
  return(t)
}

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#
# vector for chr1:22
chr<-paste0('chr', seq(1:22)) %>% sort()

chr<-'chr1'

freqTable<-
  list.files(modelDir, pattern = 'emissions_15') %>% 
  paste(modelDir, ., sep='/') %>% 
  readEmissionTable()

windowFile <- list.files(qcDir, pattern = 'dense.200.bed')
sim <- read.table(paste(qcDir, windowFile[1], sep = '/'))

sim1<-lapply(chr, function(i){
  dplyr::filter(sim, sim$V1==i) %>% 
    {freqTable[match(.[,2], rownames(freqTable)),1:(ncol(freqTable))]} %>%
    na.omit() #%>%
    #rbind(c(type, i, rep('', ncol(.)-2)), colnames(.), .)#%>%
#    write.table(paste0(simulateDir,"/", type, "_", i,"_binary.txt"), sep='\t', col.names = FALSE ,row.names = FALSE, quote = FALSE)
})

#----------------------------------------------------------------------#
# SIMULATE NOISE =======================================================
#----------------------------------------------------------------------#

set.seed(10)

## varying degrees of noise
pc<-c(1, 5, 20)
intproject<-paste0('simulate_', pc/2)

pcRemove <- apply(sim1[[1]], 2, function(x) {
  sum(x) * 0.01
}) %>% round()

index <- which(sim1[[1]][, 1] == 1) %>% unname() %>%
  sample(pcRemove[1])

sim1[[1]][index, 1] <- 0


