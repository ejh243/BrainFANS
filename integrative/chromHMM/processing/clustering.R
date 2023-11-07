##---------------------------------------------------------------------#
##
## Title: clustree =================================================
##
## Purpose of script:
##
## Author: 
##
## Date Created: 2022-08-12
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 
## set working directory
setwd("/lustre/home/jms260/BrainFANS/")

intproject<-'chromVal5_01'


source("integrative/chromHMM/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#

library(clustree)


#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#
states<-read.table(paste0(modelDir, "/robust/stateAssignments.bed")) 

#ggplot(states, aes(x=V1, y=V2)) + geom_point(stat="sum", aes(size=..n..))

pdf(file = paste0(modelDir, "/robust/cluster.pdf"),   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 8)
clustree(states, prefix = "V", node_colour = "grey")
dev.off()
#----------------------------------------------------------------------#
# CALCULATE STABILITY ==================================================
#----------------------------------------------------------------------#
# calculate stability score for each state in each model, according to proportion of
# pre-state sites going to state, and proportion of state sites coming from pre-state.

stateNo<-list.files(modelDir, pattern = 'emissions.+\\.txt') %>%
  gsub('emissions_', '', .) %>%
  gsub('.txt', '', .) %>%
  as.numeric()


stateNo<- c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32)


stateNo
stability<-NULL
for (i in 2:length(stateNo)){
  #measure is a matrix with number of states of model y 
  measure<-matrix(data=NA, ncol = stateNo[i], nrow = stateNo[i-1])
  for (x in unique(states[,i-1])){ #for each state within the model
    xtotal<-sum(states[,i-1]==x) # count number of sites in each x state 
    for (y in unique(states[,i])){
      ytotal<-sum(states[,i]==y) #count number of sites in each y state
      both<-nrow(states[states[,i-1]==(x),] %>% #calculate number of y sites of each y state are in each x state. e.g. x state==9, y state==2 88206 sites crossover
                   .[.[,i]==(y),])
      measure[x, y]<-(both/xtotal)#+(both/ytotal))/2    #scale to the amount of sites in each x or y state
    }
  }
  # scale so sum of each col (equal to all in state y is 1 
  measure<-apply(measure,2,function(x){x/sum(x)})
  stability[[i]]<-apply(measure,1, function(x){max(x)}) #save the max value (i.e. indicating the most highly similar state in the next model)  )
  ## saved in the order of the states stability[[2]][1]==state 1 in the second smallest number of state model
  
}

lapply(stability, function(x){mean(x)})

i=2
measure<-matrix(data=NA, ncol = stateNo[i], nrow = stateNo[i-1])
for (x in unique(states[,i-1])){ #for each state within the model
  xtotal<-sum(states[,i-1]==(x)) # count number of sites in each x state 
  for (y in unique(states[,i])){
    ytotal<-sum(states[,i]==(y)) #count number of sites in each y state
    both<-nrow(states[states[,i-1]==(x),] %>% #calculate number of y sites of each y state are in each x state. e.g. x state==9, y state==2 88206 sites crossover
             .[.[,i]==(y),])
    measure[x, y]<-(both/xtotal)*(both/ytotal)
  }
}
# scale so sum of each col (equal to all in state y is 1 
measure<-apply(measure,2,function(x){x/sum(x)})
stability[[i]]<-apply(measure,2, function(x){max(x)}) #save the max value (i.e. indicating the most highly similar state in the next model)  )
## saved in the order of the states stability[[2]][1]==state 1 in the second smallest number of state model



#----------------------------------------------------------------------#
# CALCULATE REPRODUCBILITY =============================================
#----------------------------------------------------------------------#


measure<-matrix(data=NA, nrow = 10, ncol = 10)
for (y in 1:10){
  xtotal<-sum(states[,1]==y) # number of 
  for (x in 1:10){
    ytotal<-sum(states[,2]==x)
    both<-nrow(states[states[,1]==(y),] %>%
                 .[.[,2]==(x),])
    measure[x, y]<-((both/xtotal)+(both/ytotal))/2
  }
}
measure<-apply(measure,2,function(x){x/sum(x)})
stability<-apply(measure,1, function(x){max(x)})
stability


x=9
corr.mat<-NULL
repro<-as.numeric(states[,1] == x)
for (y in unique(states[,2])){
  repro2<-as.numeric(states[,2] == y)
  corr.mat[x, y]<-cor.test(repro, repro2)$estimate
}


