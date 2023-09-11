##---------------------------------------------------------------------#
##
## Title: Proportion of 1s =============================================
##
## Purpose of script: to write a document with proportion of binary
##
## Author: Jessica Shields
##
## Date Created: 2022-11-23
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/home/jms260/BrainFANS/")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args[1]<-"_cmpltReadsStg1"
} 
intproject<-args[1]
project=''


source("integrative/chromHMM/config/config.r")

completeDir="//lustre/projects/Research_Project-MRC190311/ChIPSeq/sim/4_fractions/"
completeDir="/lustre/home/jms260/package/tmp/chip/subset"
completeDir="/lustre/projects/Research_Project-MRC190311/ChIPSeq/hsColorectal/4_complete/chmm/ctrl"
#completeDir="/lustre/projects/Research_Project-MRC190311/ChIPSeq/NFkB/4_complete/sample"
#completeDir = "/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/_cmpltRoadmap/"
#completeDir="/lustre/home/jms260/tmp/"


#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#

library(magrittr)
library(ggplot2)
library(gtools)
library(plyr)
library(dplyr)
library(tidyr)
library(renz)

#----------------------------------------------------------------------#
# FUNCTIONS ========================================================
#----------------------------------------------------------------------#
## adapted from the renz package.

dir.MM <- function(data, unit_S = 'mM', unit_v = 'au', plot = TRUE){
  ## ------------------ Removing incomplete data -------------------- ##
  data <- data[complete.cases(data), ]
  
  ## --------------------- Estimating the seed ---------------------- ##
  t <- ecb(data, plot = FALSE)
  K <- t$fitted_parameters[1]
  V <- t$fitted_parameters[2]
  seed = list(Km = K, Vm = V)
  
  ## ----------------------------- Fitting the curve ----------------------------- ##
  title<-colnames(data)[2]
  names(data) <- c('S', 'v')
  model <- nls(data$v ~ (Vm * data$S)/(Km + data$S), data =  data, start = seed )

  ## --------------------------- Computuing parameters --------------------------- ##
  Km <- round(summary(model)$coefficient[1,1], 3)
  sd_Km <- round(summary(model)$coefficient[1,2], 3)

  Vm <- round(summary(model)$coefficient[2,1], 3)
  sd_Vm <- summary(model)$coefficient[2,2]

  ## --------------------------- Fitted velocity lues ---------------------------- ##
  mm.eq <- function(x) {(Vm * x)/(Km + x)}

  data$fitted_v <- mm.eq(data$S)

  ## ------------------- Plotting the transformed variables ---------------------- ##
  if (plot){
    parameters <- paste(title, '     Km: ', Km, '     Vm: ', Vm, sep = "")
    plot(data$S, data$v,
         xlab = paste("[S] (", unit_S, ")", sep = ""),
         ylab = paste("v (", unit_v, ")", sep = ""),
         main = parameters, 
         ylim = c(0, Vm))
    
    x <- seq(from = 0, to = max(data$S), by = max(data$S)/1000)
    y <- mm.eq(x)
    points(x, y, ty = 'l')
    abline(v=Km, col="blue")
    abline(h= Vm, col = "red")
  }
  
  ### --------------------- Output -------------------------------- ##
  KmVm <- c(Km, Vm)
  names(KmVm) <- c("Km", "Vm")
  
  output <- list(KmVm, data)
  names(output) <- c('parameters', 'data')
  
  return(output)
}


#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#

files <- list.files(completeDir, pattern='completeness\\.', recursive = TRUE) %>%
  mixedsort()
files
#files<-files[grep('default.+', files)]
#files<-files[grep('seed', files, invert = TRUE)]
files<-files[grep('it', files, invert = TRUE)]
#files <-files[grep('_[0-3]/', files, invert = TRUE)]
#files<-files[grep('f0.02.+', files)]
files


## combine all into dataframe
comp<-NULL
for (x in seq_along(files)){
  print(x)
	tmp <- read.table(paste0(completeDir,'/', files[x]), sep='\t') %>%
		cbind(sample = c(1:nrow(.)), datatype = gsub('/.+', '',files[x]), .)
	comp <- rbind(comp, tmp)
}

comp

colnames(comp)<-c('sample', 'datatype', 'no', 'noReads', 'V1')

#comp['datatype'] <- c(rep("it_00", 6), rep ("it_01", 6), rep("it_02", 6))

## plot
pdf(paste0(completeDir, '/completenessPlot.pdf'), width = 7, height = 4)

ggplot(data=comp, aes(x=noReads, y=V1, group=datatype)) +
  geom_line(aes(color=datatype))+
  geom_point(aes(color=datatype))+
  xlab("Number of reads")+
  ylab("Proportion of marks")+
#  ylim(0, 0.16)+
  theme_bw()+
  scale_color_brewer(palette="Dark2")

dev.off()

# calculate Michaelis-Menten constants
comp$noReads <- round_any(comp$noReads, 50000000)
df <- comp %>% 
                group_by(datatype) %>% 
                filter(duplicated(datatype)|n()==1)

df <- pivot_wider(df[,c("datatype", "noReads", "V1")], names_from=datatype, values_from=V1)

## remove datatypes for which there are two values or less as not permitted in model
nonNA <- colSums(!is.na(df)) > 2
df <- df[nonNA]

pdf(paste0(completeDir, '/michaelisMenten.pdf'), width = 5, height = 4)

sapply(c(2:ncol(df)), function(x){ 
    dir.MM(df[,c(1,x)], unit_v = "proportion", unit_S = "reads") 
  })

dev.off()



#noLines <- list.files(completeDir, pattern='lines\\.', recursive = TRUE) %>%
#  mixedsort()
#sample <- NULL
#if (length(noLines) != 0){
#  for (x in seq_along(noLines)){
#    print(noLines[x])
#    tmp <- read.table(paste0(completeDir,'/', noLines[x]), sep='\t')  %>%
#      cbind(datatype = gsub('/.+', '',noLines[x]), .)
#    sample <- rbind(sample, tmp)
#  }
#}
#
#if (nrow(comp)!=nrow(sample)){
#  print('not matched')
#}
#
#comp<-cbind(comp, sample[c(1:nrow(comp)),])


## comparing between completness through chmm (replicates)
'
write.table(comp, paste0(completeDir, "/completenessAll.txt"), quote=FALSE)
ggplot(comp, aes(x=datatype, y=V1))+
  geom_boxplot()+
  theme_bw()+
  labs(title="GM15510", y = "Proportion of marks")


chmmSingle<-read.csv("/lustre/projects/Research_Project-MRC190311/ChIPSeq/NFkB/4_complete/sample/ctrl/completenessAll.csv")
ggplot(chmmSingle, aes(datatype, V1))+
  geom_bar(stat="identity")+
  theme_bw()

ggplot(chmmSingle, aes(datatype, V1))+
  geom_bar(stat="identity")+
  theme_bw()+
  labs(title="GM15510 (ctrl)", x = "Replicates", y = "Proportion of marks")
'


dir.MM(data, unit_v = "proportion", unit_S = "reads")
