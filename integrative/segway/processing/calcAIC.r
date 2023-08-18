##---------------------------------------------------------------------#
##
## Title: Calculate AIC      ========================================================
##
## Purpose of script: to calculate Akaike's Information Criterion for Segway
##
## Author: Jessica Shields
##
## Date Created: 2022-08-11
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 
## set working directory
setwd("")

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#

trainDir<-"/lustre/projects/Research_Project-MRC190311/integrative/segway/chromVal/2_train"

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#

## SEGWAY
files<-Sys.glob(file.path(paste0(trainDir, '/**/likelihood/likelihood.ll')))

likeli<-NULL
for (each in files){
  tmp<-read.table(each)
  colnames(tmp)<-('likelihood')
  tmp$states<-as.numeric(substring(each, 81, 82))
  likeli<-rbind(likeli, tmp)
}
likeli

likeli$logLikelihood<-log(likeli$likelihood)
likeli
likeli$params<-5+(likeli$states*5)+5+1+(likeli$states*5)

likeli$aic<-(2*likeli$params)-(2*likeli$likelihood)

plot(likeli$states, likeli$aic)


## CHROMHMM
logDir<-"/lustre/home/jms260/BrainFANS/integrative/chromHMM/logFiles/jms260"
files<-Sys.glob(file.path(paste0(logDir, '/**/likelihood.*.txt')))
likeli<-NULL
for (each in files){
  tmp<-read.table(each)
  colnames(tmp)<-('logLikelihood')
  tmp$states<-substring(each, 87, str_length(each)) %>%
    gsub(".txt", "", .) %>%
    as.numeric()
  likeli<-rbind(likeli, tmp)
}

likeli$params<-(likeli$states*5)+(likeli$states^2) # emission probs + transition probs
likeli$aic<-(2*likeli$params)-(2*likeli$logLikelihood)

par(mfrow=c(1,2))
plot(likeli$states, likeli$aic)

