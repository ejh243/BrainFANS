##---------------------------------------------------------------------#
##
## Title: are these two things the same?    ==========================
##
## Purpose of script:
##
## Author: 
##
## Date Created: 2022-09-30
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/home/jms260/BrainFANS/")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args<-"chromVal"
} 
intproject<-args

source("integrative/chromHMM/config/config.r")




#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#
library(ComplexHeatmap)
library(RcppHungarian)


readEmissionTable<-function(tablePath){
  t<-read.table(tablePath, sep='\t', header = TRUE, row.names = 1) %>% 
    as.matrix() %>%
    #.[, colOrder]# %>% # reorder to same col order
    {(.>=0.5)*1} # binarise dataframe (threshold 0.5)
  return(t)
}


library(lsa)



#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#
colOrder<-c("h3k27me3", "X5mc", "h3k27ac", "atac", "X5hmc")   

files<-Sys.glob(file.path(paste0(chromDir, '/**/3_model/emissions_10.txt'))) %>%
  sort()


dataDir<-"/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/"
suffix<-"/3_model/emissions_9.txt"

metric=NULL
states<-c(10)#, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
corr.all<-matrix(NA, nrow=5, ncol=11)

for (j in 1:length(states)){
  #print(j)
  suffix<-paste0("/3_model/emissions_", states[j], ".txt")
  k=1
  for (i in seq(0, 9, 2)){
    models<-c(paste0(dataDir, i, suffix), paste0(dataDir, i+1, suffix))
    # read in the two matching csvs
    for (each in models){
      emis_2<-readEmissionTable(each) 
      
      ## change to avoid NAs in final correlation matrix
      x<-which(rowSums(emis_2) == 0)
      emis_2[x, ] <- c(0, 3, 4, 5,6) # c(1, 2, 3, 4, 5) ##what ??? (when 1,2,3,4,5 it doesnt work as a true 1)
      
      if(each == models[1]){
        emis_1<-emis_2 # keep the matrix from overwriting itself
      }
    }
    
    
    corr.two<-matrix(data =NA, nrow=nrow(emis_1), ncol=nrow(emis_2))
    for (x in 1:nrow(emis_1)){
      for (y in 1:nrow(emis_2)){
        corr.two[x, y]<-cor.test(emis_1[x,], emis_2[y,], method = 'pearson')$estimate
      }
    }
    corr.two
    print(corr.two)
    #corr.two<-corr.two^2
    corr.two_bool<-corr.two==1
    
    # take the maximum correlation
    a<-apply(corr.two, 1, function(x){max(x)})
    b<-apply(corr.two, 2, function(x){max(x)})
    
    metric[[k]]<-(1-((nrow(corr.two)+ncol(corr.two))-sum(a)-sum(b))/10)
    k<-k+1
  }
  print(metric)
  corr.all[,j]<-metric
}
corr.two

heatmap.2(corr.two)

#how are they correlated? this measures the repeatability
#----------------------------------------------------------------------#
# Using euclidean distance measure   ===================================
#----------------------------------------------------------------------#
## simply the distance between two points
## for one split

euclidean <- function(a, b) sqrt(sum((a - b)^2))
colOrder<-c("h3k27me3", "X5mc", "h3k27ac", "atac", "X5hmc")   

files<-Sys.glob(file.path(paste0(chromDir, '/**/3_model/emissions_10.txt'))) %>%
  sort()

models<-c(paste0(dataDir, 0, suffix), paste0(dataDir, 1, suffix))

metric<-NULL
k=1
for (each in bootstrap){
  print(k)
  models<-c(paste0(dataDir, each, suffix), paste0(dataDir, each+1, suffix))
  emis<-NULL
  for (x in 1:2){
    emis[[x]]<-readEmissionTable(models[x])%>%
      as.data.frame()
    tmp<-row.match(matrix(emis[[x]]), uniqCombin)
    emis[[x]]<-as.matrix(emis[[x]])
    rownames(emis[[x]])<-tmp
  }
  
  emis
  
  
  euc<-matrix(data =NA, nrow=nrow(emis[[1]]), ncol=nrow(emis[[2]]))
  for (x in 1:nrow(emis[[1]])){
    for (y in 1:nrow(emis[[2]])){
      euc[x, y]<-euclidean(emis[[1]][x,], emis[[2]][y,])
    }
  }
  euc<-euc/max(euc)
  
  rownames(euc)<-rownames(emis[[1]])
  colnames(euc)<-rownames(emis[[2]])
  
  ord<-HungarianSolver(euc)$pairs
  euc<-euc[,ord[,2]]
  
  
  heatmap.2(euc, Colv=FALSE, Rowv=FALSE)
}

    
#----------------------------------------------------------------------#
# Bootstrap euclidean distance =========================================
#----------------------------------------------------------------------#
states<-c(5)# 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
bootstrap<-(seq(0, 9, 2))#, seq(22, 40, 2))
euc.all<-matrix(NA, nrow=length(states), ncol=length(bootstrap))

## internal for loop reads in the two emission probabilities and 
for (z in 1:length(states)){
  metric<-NULL
  suffix<-paste0("/3_model/emissions_", states[z], ".txt")
  k=1
  for (each in bootstrap){
    print(k)
    models<-c(paste0(dataDir, each, suffix), paste0(dataDir, each+1, suffix))
    emis<-NULL
    #print(models)
    
    for (x in 1:2){
      emis[[x]]<-readEmissionTable(models[x])
    }
        
    euc<-matrix(data =NA, nrow=nrow(emis[[1]]), ncol=nrow(emis[[2]]))
    for (x in 1:nrow(emis[[1]])){
      for (y in 1:nrow(emis[[2]])){
        euc[x, y]<-euclidean(emis[[1]][x,], emis[[2]][y,])
      }
    }
    euc<-euc/max(euc)
    rownames(euc)<-seq(1:nrow(emis[[1]]))
    colnames(euc)<-seq(1:nrow(emis[[2]]))
    
    
    ord<-HungarianSolver(euc)$pairs
    euc<-euc[,ord[,2]]
    
    
    heatmap.2(euc, Colv=FALSE, Rowv=FALSE)
    
    metric[[k]]<-apply(euc, 2, min) %>%
      sum()
    print(metric)
    k<-k+1
  }
  euc.all[z,]<-metric
}
## matrix where the cols are bootstrap and rows are states
rownames(euc.all)<-paste0('E', seq(states[1],states[length(states)]))

reshape2::melt(euc.all) %>%
  ggplot(aes(Var1, value))+
  geom_violin()+
  stat_summary(fun="mean", 
               geom="point", color="black")+
  xlab("Number of States in Model")+
  theme_bw()
#----------------------------------------------------------------------#
#             ==========================================================
#----------------------------------------------------------------------#

i=0
models<-c(paste0(dataDir, i, suffix), paste0(dataDir, i+1, suffix))
# read in the two matching csvs
for (each in models){
  emis_2<-readEmissionTable(each)
  
  ## change to avoid NAs in final correlation matrix
  #x<-which(rowSums(emis_2) == 0)
  #emis_2[x, ] <- c(0, 3, 4, 5,6) # c(1, 2, 3, 4, 5) ##what ??? (when 1,2,3,4,5 it doesnt work as a true 1)
  
  if(each == models[1]){
    emis_1<-emis_2 # keep the matrix from overwriting itself
  }
}


corr.two<-matrix(data =NA, nrow=nrow(emis_1), ncol=nrow(emis_2))
for (x in 1:nrow(emis_1)){
  for (y in 1:nrow(emis_2)){
    corr.two[x, y]<-cor.test(emis_1[x,], emis_2[y,], method = 'pearson')$estimate
  }
}
corr.two<-cor()
print(corr.two)
#corr.two<-corr.two^2
corr.two_bool<-corr.two==1

# take the maximum correlation
a<-apply(corr.two, 1, function(x){max(x)})
b<-apply(corr.two, 2, function(x){max(x)})

metric[[k]]<-(1-((nrow(corr.two)+ncol(corr.two))-sum(a)-sum(b))/10)


#----------------------------------------------------------------------#
# same but no state split ==============================================
#----------------------------------------------------------------------#
## only for a 10 state model
k=1
for (i in seq(0, 9, 2)){ ## for bootstraps 0-9 (five data splits)
  #models<-c(files[i],files[i+1])
  models<-c(paste0(dataDir, i, suffix), paste0(dataDir, i+1, suffix))
  # read in the two matching csvs
  for (each in models){
    emis_2<-readEmissionTable(each) # reorder to same col order
    
    ## change to avoid NAs in final correlation matrix
    x<-which(rowSums(emis_2) == 0)
    emis_2[x, ] <- c(0, 3, 4, 5,6) # c(1, 2, 3, 4, 5) ##what ??? (when 1,2,3,4,5 it doesnt work as a true 1)
    
    if(each == models[1]){
      emis_1<-emis_2 # keep the matrix from overwriting itself
    }
  }
  
  
  corr.two<-matrix(data =NA, nrow=nrow(emis_1), ncol=nrow(emis_2))
  for (x in 1:nrow(emis_1)){
    for (y in 1:nrow(emis_2)){
      corr.two[x, y]<-cor.test(emis_1[x,], emis_2[y,], method = 'pearson')$estimate
    }
  }
  corr.two
  corr.two_bool<-corr.two==1
  
  print(corr.two)
  
  # take the maximum correlation
  a<-apply(corr.two, 1, function(x){max(x)})
  b<-apply(corr.two, 2, function(x){max(x)})
  
  metric[[k]]<-(1-((nrow(corr.two)+ncol(corr.two))-sum(a)-sum(b))/10)
  k<-k+1
}
metric
boxplot(metric)

as.data.frame(metric) %>%
  ggplot(., aes(y=metric))+
  geom_boxplot()+
  theme_bw()
## therefore on average there is a 0.93 correlation between the two splits of the data  


i=0
models<-c(paste0(dataDir, i, suffix), paste0(dataDir, i+1, suffix))
# read in the two matching csvs
for (each in models){
  emis_2<-readEmissionTable(each)
  
  ## change to avoid NAs in final correlation matrix
  x<-which(rowSums(emis_2) == 0)
  emis_2[x, ] <- c(0, 3, 4, 5,6) # c(1, 2, 3, 4, 5) ##what ??? (when 1,2,3,4,5 it doesnt work as a true 1)
  
  if(each == models[1]){
    emis_1<-emis_2 # keep the matrix from overwriting itself
  }
}

emis_1
corr.two<-matrix(data =NA, nrow=nrow(emis_1), ncol=nrow(emis_2))
for (x in 1:nrow(emis_1)){
  for (y in 1:nrow(emis_2)){
    corr.two[x, y]<-cor.test(emis_1[x,], emis_2[y,], method = 'pearson')$estimate
  }
}
corr.two
corr.two_bool<-corr.two==1

print(corr.two)

# take the maximum correlation
a<-apply(corr.two, 1, function(x){max(x)})
b<-apply(corr.two, 2, function(x){max(x)})




plot(emis_1[which(a<1), ])
plot(emis_2[which(b<1),])

#----------------------------------------------------------------------#
# plotted graphic with peaks ===========================================
#----------------------------------------------------------------------#

## plotted graphic with peaks 
ggplot(emis_2[which(b<1),1])
p<-NULL
for (x in 1:5){
  binary<-emis_2[which(b<1),x]
  if (binary==0){
    data<-data.frame(x = 1:5, y = 0)
    p[[x]]<-ggplot(data, aes(x,y))+
      geom_line()+
      ylim(0, 1)+
      theme_void()
  } else {
    data<-data.frame(x = 1:5, y = c(0, 0, 1, 0, 0))
    p[[x]]<-ggplot(data, aes(x,y))+
      geom_line()+
      theme_void()
  }
}
x=6
for (y in 1:5){
  binary<-emis_1[which(a<1),y]
  if (binary==0){
    data<-data.frame(x = 1:5, y = 0)
    p[[x]]<-ggplot(data, aes(x,y))+
      geom_line()+
      ylim(0, 1)+
      theme_void()
  } else {
    data<-data.frame(x = 1:5, y = c(0, 0, 1, 0, 0))
    p[[x]]<-ggplot(data, aes(x,y))+
      geom_line()+
      theme_void()
  }
  x<-x+1
}
grid.arrange(grobs = p,ncol=2, as.table=FALSE)

data<-data.frame(x = 1:5, y = 0)
ggplot(data, aes(x,y))+
  geom_line()+
  ylim(0, 1)
  theme_void()

data<-data.frame(x = 1:5, y = c(0, 0, 1, 0, 0))
ggplot(data, aes(x,y))+
  geom_line()+
  theme_void()

grid.arrange()

#----------------------------------------------------------------------#
# heatmap graphic of conserved =========================================
#----------------------------------------------------------------------#

## figure out where the issue is AND WHY??
for (each in files){
  tmp<-read.table(each, sep='\t', header = TRUE)
  if (ncol(tmp)<6){
    print(each)
    print(dim(tmp))
  }
}


#----------------------------------------------------------------------#
#   distribution of segment lengths   ==================================
#----------------------------------------------------------------------#
## the distribution of segment lengths
## mean(log(x)) rather than log(mean(x)) https://stats.stackexchange.com/questions/250209/log-mean-vs-mean-log-in-statistics


dt<-fread(paste0(modelDir,"/N+_10_segments.bed")) %>%
  cbind(.[,3]-.[,2]) %>%
  cbind(log(.[,5]))

dt<-dt[,c(4,5,6)] %>%
  as.data.frame() %>%
  setnames(c('state', 'length', 'log'))




df<-dt %>% dplyr::group_by(state) %>% 
  dplyr::summarise(
#  average = mean(length),
#  standDev = sd(length),
    logAverage = mean(log),
    logSD = sd(log),
    logMed = median(log)
)

df$state<-mixedsort(gsub('E', '', df$state))

df$state<-rownames(emis[[nrow(df)-4]])

ggplot(df, aes(state, logMed))+
  geom_point()+
#  geom_errorbar(aes(ymin=df$logMed-df$standDev/2, ymax=df$average+df$standDev/2))+
  theme_bw()

ggplot(df, aes(state, logAverage))+
  geom_point()+
  #  geom_errorbar(aes(ymin=df$average-df$standDev/2, ymax=df$average+df$standDev/2))+
  theme_bw()



ggplot(df, aes(state, logAverage))+
  geom_boxplot()+
  theme_bw()
  stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2)

df2<-dt %>% dplyr::group_by(state) %>% 
  dplyr::summarise(
    #  average = mean(length),
    #  standDev = sd(length),
    logAverage = mean(log),
    logSD = sd(log),
    logMed = median(log)
  )

df2$state<-mixedsort(gsub('E', '', df2$state))

df2$state<-rownames(emis[[nrow(df2)-4]])

df<-rbind(df, df2)

###################
x=5
dt<-fread(paste0(modelDir,"/N+_",x,"_segments.bed")) %>%
  cbind(.[,3]-.[,2]) %>%
  cbind(log(.[,5]))

dt<-dt[,c(4,5,6)] %>%
  as.data.frame() %>%
  setnames(c('state', 'length', 'log'))

dt$newstate<-gsub('E', '', dt$state)
dt$newstate<-factor(dt$newstate, levels =  sort(rownames(emis[[x-4]]), decreasing = TRUE))

plots[[4]]<-ggplot(dt, aes(length,newstate))+
  geom_violin()+
  theme_bw()+
  xlab('')+
  ylab('')+
  coord_trans(x="log10")
  
grid.arrange(grobs = plots, ncol=2)

x=6
dt<-fread(paste0(modelDir,"/N+_",x,"_segments.bed")) %>%
  cbind(.[,3]-.[,2]) #%>%
#  cbind(log(.[,5]))

dt<-dt[,c(4,5)] %>%
  as.data.frame() %>%
  setnames(c('state', 'length'))#, 'log'))

dt$state<-gsub('E', '', dt$state)


ord<-rownames(emis[[x-4]])
for (y in 1:length(ord)){
  print(as.character(y))
  dt$state[dt$newstate == as.character(y)] <- ord[y]
}
head(dt)

dt$newstate<-factor(dt$state, levels =  order(data.frame(rownames(emis[[x-4]]))[,1]))

plots[[1]]<-ggplot(dt, aes(newstate, length))+
  geom_violin()+
  theme_bw()+
  xlab('')+
  ylab('')+
#  scale_y_continuous(trans='log2')+
#  scale_y_reverse()
#  scale_y_log10()
  coord_trans(y="log10")

plots[[2]]<-ggplot(NULL)+theme_void()
grid.arrange(grobs = plots, ncol=2)
####################

## addition for split02
df2<-fread("/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/1/3_model/N+_10_segments.bed")
df2$V5<-df2$V3-df2$V2
df2$V6<-log(df2$V5)

ggplot(df2, aes(V4, V6))+
  geom_boxplot()+
  theme_bw()

ggplot(df2, aes(V4, V6))+
  geom_violin()+
  #geom_boxplot(width=0.05)+
  theme_bw()


plot_grid(plotlist = plots, align='hv', axis = "lrbt", greedy = TRUE)
#----------------------------------------------------------------------#
# The genomic coverage of each label  ==================================
#----------------------------------------------------------------------#
cov<-read.table("/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/0/3_model/N+_10_overlap.txt", sep='\t', header=TRUE)
cov<-cov[c(1:10),c(1,2)]
colnames(cov)<-c('State', 'split01')
split02<-read.table("/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/bootstrap/1/3_model/N+_10_overlap.txt", sep='\t', header=TRUE) %>%
  .[c(1:10),2]


## for the increasing number of states on the whole genome
set.seed(150)
brbg <- c(brewer.pal(12, 'Set3'), brewer.pal(3, 'Dark2')) %>%
  sample()

cov<-matrix(NA, nrow=110, ncol = 3) #110, 12
State<-data.frame(rownames(emis[[1]]))
State<-rownames(emis[[1]])
for (x in 2:length(emis)){
  State<-c(State,rownames(emis[[x]]))
}
cov[,1]<-State

# read in each of the genomic coverages and add to long_pivoted df
x=5
cov[1:x,2]<-paste0('E', x)
cov[1:x,3]<-read.table(paste0(modelDir, "/N+_",x,"_overlap.txt"), sep='\t', header=TRUE) %>%
  .[c(1:x),2]
y<-1+x
for (x in 6:15){
  cov[y:(y+x-1),2]<-paste0('E', x)
  cov[y:(y+x-1),3]<-read.table(paste0(modelDir, "/N+_",x,"_overlap.txt"), sep='\t', header=TRUE) %>%
    .[c(1:x),2] 
  y<-y+x
}
#colnames(cov)<-c('State', 'variable', 'value')
cov

cov<-as.data.frame(cov) %>%
  setnames(c('State', 'Model', 'genomeCov'))
cov$genomeCov<-as.numeric(levels(cov$genomeCov))[cov$genomeCov]

cov[,2]<-factor(cov[,2], levels = unique(cov[,2]))

set.seed(120)
brbg <- c(brewer.pal(9,"YlOrRd"), brewer.pal(2, 'Dark2')) %>%
  sample()


ggplot(cov, aes(fill = State, x = Model, y = genomeCov))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.7, colour = "black") +
  scale_fill_manual(values = brbg)+
  theme_bw()
## the enrichment of each label for previously annotated genomic elements


#----------------------------------------------------------------------#
#             ==========================================================
#----------------------------------------------------------------------#

tmp<-df[which(df[,4]!="E4"),]

ggplot(tmp, aes(V4, V5))+
  geom_violin()+
  ylim(0,1000)+
  theme_bw()

tmp<-df[which(df[,4]=="E4"),]

ggplot(tmp, aes(V4, V5))+
  geom_hist()

hist(tmp$V6)

#----------------------------------------------------------------------#
# Comparing matchingness bn the two splits of the dataset across states 
#----------------------------------------------------------------------#
states<-c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
bootstrap<-c(seq(0, 18, 2))#, seq(22, 40, 2))
corr.all<-matrix(NA, nrow=length(states), ncol=length(bootstrap))

## internal for loop reads in the two emission probabilities and 
for (z in 1:length(states)){
  metric<-NULL
  suffix<-paste0("/3_model/emissions_", states[z], ".txt")
  y=1
  for (each in bootstrap){
    models<-c(paste0(dataDir, each, suffix), paste0(dataDir, each+1, suffix))
    emis<-NULL
    
    for (x in 1:2){
      emis[[x]]<-readEmissionTable(models[x]) %>%
        split(seq(nrow(.))) # convert matrix rows into single list entry
    }
    
    compare<-cbind(emis[[1]] %in% emis[[2]], 
                   emis[[2]] %in% emis[[1]])
    metric[[y]]<-(1-(length(compare)-sum(compare[,1])-sum(compare[,2]))/10)
    y<-y+1
  }
  corr.all[z,]<-metric
}
## matrix where the cols are bootstrap and rows are states
rownames(corr.all)<-paste0('E', seq(states[1],states[length(states)]))


## PLOT as bixplot
reshape2::melt(corr.all) %>%
  ggplot(aes(Var1, value))+
  geom_violin()+
  stat_summary(fun="mean", 
               geom="point", color="black")+
  xlab("Number of States in Model")+
  theme_bw()


## plot points with error bars instead of boxplots
mean<-apply(corr.all, 1, function(x){mean(x)})
sd<-apply(corr.all[,1:10], 1, function(x){sd(x)})
corr.all<-cbind(corr.all, mean, sd)
ggplot(corr.all, aes(rownames(corr.all), mean))+
  geom_point()+
  xlab("Number of States in Model")+
  theme_bw()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))

#----------------------------------------------------------------------#
# Comparing correlation across all  ====================================
#----------------------------------------------------------------------#
states<-c(5, 6) #, 7, 8, 9, 10, 11, 12, 13, 14, 15)
plots<-NULL
noState<-matrix(nrow=length(states), ncol=2)
noState[,1]<-states
for (z in 1:length(states)){
  suffix<-paste0("/3_model/emissions_", states[z], ".txt")
  emis<-NULL
  filePath<-paste0(dataDir, 0, suffix)
  emis[[1]]<-readEmissionTable(filePath) %>%
    split(seq(nrow(.))) # convert matrix rows into single list entry
  for (x in c(1:19,22:41)){
    filePath<-paste0(dataDir, x, suffix)
    #print(filePath)
    emis[[x]]<-readEmissionTable(filePath) %>%
      split(seq(nrow(.))) # convert matrix rows into single list entry
    emis[[1]]<-Reduce(intersect,list(emis[[1]],emis[[x]]))
    }
  #print(length(emis[[1]]))#
  noState[z,2]<-length(emis[[1]])
  consvd<-emis[[1]] %>%
    as.data.frame(., col.names=paste0("E", seq(len=length(.))), row.names = colOrder)%>%
    t()
  
  print(Heatmap(consvd, col = brewer.pal(10,"YlOrRd"), show_row_dend = FALSE, 
                      show_column_dend =FALSE,show_heatmap_legend = FALSE))
}
noState
emis

eval(plots$call)


consvd<-emis[[1]] %>%
  as.data.frame(., col.names=paste0("E", seq(len=length(.))), row.names = colOrder)%>%
  t()

heatmap.2(consvd, scale = "none", Rowv=FALSE, Colv = FALSE, 
             trace="none", col = brewer.pal(10,"YlOrRd"), keysize = 0.1, key=FALSE, 
             margins = c(9,9))

Heatmap(consvd, col = brewer.pal(10,"YlOrRd"), show_row_dend = FALSE, show_column_dend =FALSE,show_heatmap_legend = FALSE)


p<-Heatmap(consvd)

p+Heatmap(consvd)

##but do i actually want to have a proportion assigned to the likelihood of finding that state, it may not show up every single time but 
plot(noState)

#----------------------------------------------------------------------#
# Generate a matrix of all combinations, then ask how often each unique
# combination is found in the grid. 
#----------------------------------------------------------------------#
states<-c(5) #, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
noState<-matrix(nrow=length(states), ncol=2)
noState[,1]<-states
#  suffix<-paste0("/3_model/emissions_", states[z], ".txt")
emis<-NULL
filePath<-paste0(dataDir, 0, suffix)
emis<-readEmissionTable(filePath) %>%
  split(seq(nrow(.)))%>%
  sapply(paste0, collapse="")%>%
  as.data.frame()
for (x in c(1:19,22:41)){
  filePath<-paste0(dataDir, x, suffix)
  emis<-rbind(emis, readEmissionTable(filePath) %>%
    split(seq(nrow(.)))%>%
    sapply(paste0, collapse="")%>%
    as.data.frame())
}

emis
uniqCombin<-unique(emis)
for(x in 1:nrow(uniqCombin)){
  uniqCombin[x,2]<-length(which(emis[,1]==uniqCombin[x,1]))
}
uniqCombin

for(x in 1:nrow(uniqCombin)){
  uniqCombin[x,2]<-length(which(emis[,1]==uniqCombin[x,1]))
}

uniqCombin[,2]/40
sum(uniqCombin$V2)

#----------------------------------------------------------------------#
#    get unique combinations per state???===============================
#----------------------------------------------------------------------#
##get unique combinations per state???

states<-c(5) #, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
noState<-matrix(nrow=length(states), ncol=2)
noState[,1]<-states
#  suffix<-paste0("/3_model/emissions_", states[z], ".txt")
emis<-NULL
filePath<-paste0(dataDir, 0, suffix)
emis<-readEmissionTable(filePath) %>%
  as.list()
for (x in c(1:19,22:41)){
  filePath<-paste0(dataDir, x, suffix)
  emis<-rbind(emis, readEmissionTable(filePath) %>%
 #               split(seq(nrow(.)))%>%
#                sapply(paste0, collapse="")%>%
                as.data.frame())
}

uniqCombin<-unique(emis)
uniqCombin
rownames(uniqCombin)<-seq(1, nrow(uniqCombin))


#----------------------------------------------------------------------#
#             ==========================================================
#----------------------------------------------------------------------#
emis[[1]]<-readEmissionTable(filePath) %>%
  split(seq(nrow(.)))%>%
  sapply(paste0, collapse="")%>%
  as.data.frame()


#., col.names=paste0("E", seq(len=length(.))), row.names = colOrder)



print(length(emis[[1]]))#
noState[z,2]<-length(emis[[1]])

consvd<-Reduce(intersect,list(emis[[1]],emis[[2]])) %>%
  as.data.frame(., col.names=paste0("E", seq(len=length(.))), row.names = colOrder)%>%
  t()

heatmap.2(consvd, scale = "none", Rowv=FALSE, Colv = FALSE, 
          trace="none", col = brewer.pal(10,"YlOrRd"), keysize = 0.1, key=FALSE, 
          margins = c(9,9))


split1<-split(emis_1, seq(nrow(emis_1)))
split2<-split(emis_2, seq(nrow(emis_2)))
corr.1<-NULL
corr.2<-NULL
for (x in 1:length(split2)){
    corr.1[[x]]<-split1[x] %in% split2
    corr.2[[x]]<-split2[x] %in% split1
}
data<-cbind(as.numeric(corr.1), as.numeric(corr.2))

metric1<-(1-(length(data)-sum(data[,1])-sum(data[,2]))/10)
metric[[k]]<-metric1

## conserved mark sets
consvd<-Reduce(intersect,list(split1,split2)) %>%
  as.data.frame(., col.names=paste0("E", seq(len=length(.))), row.names = colOrder)%>%
  t()

heatmap.2(consvd, scale = "none", Rowv=FALSE, Colv = FALSE, 
          trace="none", col = brewer.pal(10,"YlOrRd"), keysize = 0.1, key=FALSE, 
          margins = c(9,9))


#----------------------------------------------------------------------#
#             ==========================================================
#----------------------------------------------------------------------#
emis<-NULL
suffix<-paste0("/3_model/emissions_5.txt")
filePath<-paste0(dataDir, 0, suffix)
emis<-readEmissionTable(filePath) %>%
  split(seq(nrow(.)))%>%
  sapply(paste0, collapse="")%>%
  as.data.frame()
for (x in 1){
  filePath<-paste0(dataDir, x, suffix)
  emis<-cbind(emis, readEmissionTable(filePath) %>%
                split(seq(nrow(.)))%>%
                sapply(paste0, collapse="")%>%
                as.data.frame())
}
emis<-apply(emis, 2, as.numeric)

emis<-cbind(emis, emis[match(emis[,1], emis[,2]),2])


for (x in 1:nrow(emis)){
  print(emis[x,1])
  for (y in 1:nrow(emis)){
    print(emis[y,2])
    corr.two[x, y]<-cor.test(emis[x,1], emis_2[y,2], method = 'pearson')$estimate
  }
}
corr.two

#----------------------------------------------------------------------#
# TRANSITION PROBABILITY   =============================================
#----------------------------------------------------------------------#
Tfiles<-list.files(modelDir, pattern = 'transitions_..?.txt') %>%
  mixedsort()
Efiles<-list.files(modelDir, pattern = 'emissions_..?.txt') %>%
  mixedsort()
  
  
emis<-NULL
trans<-NULL
for (x in 1:length(Efiles)){
  emis[[x]]<-readEmissionTable(paste(modelDir, Efiles[x], sep = '/'))%>%
    as.data.frame()
  
  # create list of states matched to the unique combinations
  tmp<-row.match(matrix(emis[[x]]), uniq)
  
#  emis[[x]]<-as.matrix(emis[[x]])
  #  rownames(emis[[x]])<-tmp
  
  trans[[x]]<-read.table(paste(modelDir, Tfiles[x], sep='/'), sep='\t', header = TRUE, row.names = 1) %>% 
    as.matrix()
  
  rownames(trans[[x]])<-tmp
  colnames(trans[[x]])<-tmp
}


transProb<-NULL
trans2<-NULL
for (each in trans){
  print(sum(each))
  tmp<-NULL
  for (x in 1:nrow(each)){
    tmp[[x]]<-sum(each[,x])-each[x,x]
  }
  transProb[[x]]<-(sum(tmp)/nrow(each))
  trans2[[x]]<-tmp
  names(trans2[[x]])<-rownames(each)
#  heatmap.2(each, Rowv = FALSE, Colv = FALSE, col = viridis, keysize = 0.1, key = FALSE)
}

plot(transProb)

df<-data.frame(unlist(trans2))
df$state<-names(unlist(trans2))
df

ggplot(df, aes(state, unlist.trans2.))+
  geom_violin()

tmp<-NULL
for (x in 1:nrow(each)){
  tmp[[x]]<-sum(each[,x])-each[x,x]
}
print(sum(tmp)/nrow(each))

lapply(trans, function(x){
  sum(trans[,1])-trans[1,1]
})