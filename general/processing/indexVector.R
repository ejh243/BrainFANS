##---------------------------------------------------------------------#
##
## Title: Average lines ================================================
##
## Purpose of script: to calculate average number of reads in dataset
##
## Author: Jessica Shields
##
## Date Created: 2022-11-08
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list = ls())
## set working directory


args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <-c('ChIPSeq/epiGaba', "H.276.GLU.27ac", "H.276.SOX.27ac", "H.286.GLU.27ac", "H.286.SOX.27ac", "H.344.SOX.27ac", "H.372.GLU.27ac", "H.372.SOX.27ac")
}

project<-args[1]
intproject<-''
files <- args[2:length(args)]

source("/lustre/home/jms260/BrainFANS/integrative/chromHMM/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#

library(plyr)
library(magrittr)

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#
alignQC<-
  read.table(
    paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt"), 
    sep = "\t", 
    header = TRUE, 
    stringsAsFactors = FALSE)

matched<-NULL
for (each in files){
  matched<-sapply(alignQC$Sample, function(x){match(substr(x, 1, nchar(each)), each)}) %>%
        .[!is.na(.)] %>%
        c(.,matched)
}

avLines <- alignQC[which(alignQC$Sample %in% names(matched)),"paired_total"] %>% #filter to input file rows and only total seq col
    mean() %>%
    round_any(1000000)

options("scipen"=10, "digits"=7)
cat(avLines)
