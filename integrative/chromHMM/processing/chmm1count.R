##---------------------------------------------------------------------#
##
## Title: Proportion of 1s =============================================
##
## Purpose of script: to write a document with proportion of binary
##
## Author: Jessica Shields
##
## Date Created: 2022-11-10
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/home/jms260/BrainFANS/")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args[1]<-"ChIPSeq/NFkB"
  args[2]<-"ctrl/GM15510_0"
} 
project<-args[1]
mark<-args[2]

source("integrative/chromTools/config/config.r")

#dir.create(paste0(outDir,'/', args[3]), recursive=TRUE)

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#

files <- list.files(chromToolDir, pattern = 'binary.+')
allChr <- NULL
for (each in files) {
  tmp <- read.delim(paste(chromToolDir, each, sep = '/'), skip = 1)
  allChr <- rbind(allChr, tmp)
}
prop <- apply(allChr, 2, function(x){
  sum(x)/nrow(allChr)
  })

write.table(prop[1], file = paste0(chromToolDir, "/completeness.txt"), 
            sep = '\t', col.names = FALSE, 
            row.names = FALSE)


