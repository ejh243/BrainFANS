##---------------------------------------------------------------------#
##
## Title: Identify optimal number of chromatin states
##
## Purpose of script: to identify the minimal number of states in a ChromHMM model that capture
##                    the combinations of marks present in the data
##
##                    Strategies derived from: https://www.nature.com/articles/s41586-020-2093-3
##
## Author: Jessica Shields
##
## Date Created: 2022-03-04
##
##---------------------------------------------------------------------#
## clear the R environment
rm(list=ls()) 

## set working directory
modelDir<-"/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/chromMVal/4_model/"
  setwd("/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/chromMVal/4_model/")

col<- wes_palette("Darjeeling1", 15, type = "continuous")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(tidyverse)
library(magrittr)
library(ggfortify)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
corrModel<- read.table('compare/comparedModel.txt', header = TRUE, sep = '\t')
## the maximum correlation of each state in the selected model between its best matching state in each other model

# remove reference model  
corrModel<- corrModel[, colSums(corrModel) != nrow(corrModel)]
corrModel

#the emission probabilities from all the 23 models
emisFiles <- list.files(modelDir, pattern = "emissions.*\\.txt")
emisALL<-NULL
model<-NULL
state<- NULL
for(each in emisFiles){
  model<-c(model, (str_sub(each, 11, -5)))
  tmp<-read.table(paste0(modelDir, each), header=TRUE, sep='\t', row.names = 1)
  state<-c(state, rownames(tmp))
  rownames(tmp)<-paste0((str_sub(each, 11, -5)), '_', rownames(tmp))
  emisALL<-rbind(emisALL, tmp)
}
#emisALL <- t(emisALL)
emisALL

# sort vector containing the state model numbers
model<-sort(as.numeric(model))

emisMETA<- as.data.frame(state, row.names = rownames(emisALL))
#----------------------------------------------------------------------#
# Strategy 1
#----------------------------------------------------------------------#
## 1. calculate the median correlation of all the states in the reference model against the simpler models
corrModel$rowMedian = apply(corrModel[,-1], 1, median)

## 2. plot these numbers against the number of states in the model
ggplot(data=corrModel, aes(x=State, y=rowMedian)) +
  geom_point()+
  theme_bw()+
  scale_x_continuous("State", labels = as.character(corrModel$State), breaks = corrModel$State) 


# 3. look at the number of states at which both series reached a plateau


#----------------------------------------------------------------------#
# Strategy 2
#----------------------------------------------------------------------#
## As a complementary strategy, the emission probabilities from all the 23 
## models (considering 2–24 states) from both replicates were clustered together. 
## The rationale behind this strategy is that very similar states across models 
## will tend to cluster together, so there must be an optimal number of clusters 
## corresponding to the optimal number of states in the model
emisALL<-scale(emisALL)

## k-means clustering with k between 2 and 24
all<-NULL
for (x in 1:100){
  tmp<-NULL
  for (y in min(model):max(model)){
    kmean<- kmeans(emisALL, y)
    tmp<-rbind(tmp, kmean$betweenss/kmean$totss) 
  }
  all<-cbind(all,tmp)
}

all<- as.data.frame(all)
all$model<- model
all$average<- apply(all, 1, mean)
all$ratio<- all[,'average']/all[nrow(all), 'average']  

ggplot(data=all, aes(x=model, y=ratio)) +
  geom_point()+
  geom_abline(intercept = 0.95, slope = 0, linetype="dashed")+
  theme_bw()+
  scale_x_continuous("Model", labels = as.character(corrModel$State), breaks = corrModel$State)

kmean<- kmeans(emisALL, 11)
autoplot(kmean, emisALL, frame=TRUE)+
  theme_bw()

