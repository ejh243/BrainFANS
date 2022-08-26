##---------------------------------------------------------------------#
##
## Title: Make Group Analysis File ====================================
##
## Purpose of script: to create a txt metadata file with the n+ n- and t for cell fractions
##
## Author: Jessica Shields
##
## Date Created: 2022-07-27
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/projects/Research_Project-MRC190311/scripts/sequencing")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args[1]<-"DNAhydroxy/MRC"
  args[1]<-'WGBS/rizzardi'
  args[1]<-'ATACSeq/MRC'
  #args[1]<-'ChIPSeq/epiGaba'
  args[2]<-'prefrontal cortex|PFC'
} 

project<-args[1]
tissue<-paste(args[2:length(args)], collapse = ' ')

print(tissue)

source("BSSeq/config/config.r")

#----------------------------------------------------------------------#
# LOAD DATA ==========================================================
#----------------------------------------------------------------------#

library(magrittr)
library(rlang)

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#
sampleSheet<-read.csv(sampleSheet)

if (file.exists(paste0(metaDir, '/stage1Samples.txt'))){
  sampleid<-read.table(paste0(metaDir, '/stage1Samples.txt'))[,1]
} else {
  print('Using unfiltered samples.txt')
  sampleid<-read.table(paste0(metaDir, '/samples.txt'))[,1]
}

sampleSheet<-sampleSheet[which(sampleSheet$sampleID %in% sampleid),]

## filter sampleSheet by selected tissue
if (length(sampleSheet$tissue)!=0){
  sampleSheet<-sampleSheet[grep(tissue, sampleSheet$tissue), ]
}


# create reference dictionary to rename fractions
#dic<-data.frame( 
#  c('N-', 'olig|sox10|neun neg'),
#  c('N+', 'glu|gaba|neun|neun pos'),
#  c('T', 'total|bulk'), 
#  c('DN', 'doubleneg')
#)

# rename fractions to N-, N- and S
cell <- toupper(sampleSheet$fraction)
cell
#for (x in 1:length(colnames(dic))){
#  cell<-replace(cell, 
#                grepl(dic[2,x], cell), 
#                as.character(dic[1,x]))
#}


## get assay target
mark<-as.character(sampleSheet$target)

mark

samples<-cbind(as.character(sampleSheet$sampleID), cell, mark)

if (is_empty(sampleSheet$controlID) ==FALSE){
  samples<-cbind(samples, as.character(sampleSheet$controlID))
}

write.table(samples, paste0(metaDir, '/samplesForGroupAnalysis.txt'),
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
