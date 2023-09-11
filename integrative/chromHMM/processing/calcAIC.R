##---------------------------------------------------------------------#
##
## Title: Calculate AIC      ========================================================
##
## Purpose of script: to calculate Akaike's Information Criterion for ChromHMM models
##
## Author: Jessica Shields
##
## Date Created: 2022-08-25
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/projects/Research_Project-MRC190311/scripts/")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args<-"simulate_1"
} 
intproject<-args

source("integrative/chromHMM/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#

library(ggplot2)
library(magrittr)

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#

files<-Sys.glob(file.path(paste0(qcDir, '/likelihood.*.txt')))

likeli<-NULL
for (each in files){
  if (file.info(each)$size != 0){ #if file not empty
    tmp<-read.table(each)
    colnames(tmp)<-('logLikelihood')
    tmp$states<-gsub(paste0(qcDir, '/likelihood.'), '',each) %>%
        gsub(".txt", "", .) %>%
        as.numeric()
    likeli<-rbind(likeli, tmp)
  }
}

likeli$params<-(likeli$states*5)+(likeli$states^2) # emission probs + transition probs
likeli$aic<-(2*likeli$params)-(2*likeli$logLikelihood)

ggplot(likeli, aes(states, aic))+
  geom_point()+
  theme_bw()+
  geom_line()


likeli

write_csv(likeli, paste0(modelDir, "/robust/aic.csv"))


# save plot to png
png(paste0(modelDir, "/robust/aic.png"))
ggplot(likeli, aes(states, aic))+
  geom_point()+
  theme_bw()+
  geom_line()
dev.off()
