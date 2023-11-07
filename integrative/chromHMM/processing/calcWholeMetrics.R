##---------------------------------------------------------------------#
##
## Title: State quality metrics for complete dataset =========
##
## Purpose of script:
##
## Author: 
##
## Date Created: 2022-10-12
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/home/jms260/BrainFANS/")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args<-"chromValAll"
} 
intproject<-args

source("integrative/chromHMM/config/config.r")

colOrder<-c("h3k27me3", "h3k27ac", "atac", "X5hmc")

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#

library(gtools)
library(gridGraphics)
library(grid)
library(tidyfast)

readEmissionTable<-function(tablePath){
  t<-read.table(tablePath, sep='\t', header = TRUE, row.names = 1) %>% 
    as.matrix() %>%
    .[, colOrder] %>% # reorder to same col order
  {(.>=0.5)*1} # binarise dataframe (threshold 0.5)
  return(t)
}

euclidean <- function(a, b) sqrt(sum((a - b)^2))

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#

files<-list.files(modelDir, pattern = 'emissions_..?.txt') %>%
  mixedsort()
files

#----------------------------------------------------------------------#
# 1. Generate vector with unique combinations of marks across all models 
#----------------------------------------------------------------------#
uniq<-NULL
for (each in files){
  uniq<-rbind(uniq, readEmissionTable(paste(modelDir, each, sep = '/')) %>%
                as.data.frame())
}

uniq<-unique(uniq)
rownames(uniq)<-seq(1, nrow(uniq))
uniq

#----------------------------------------------------------------------#
# 2. Heatmaps of combinations of marks as number states increases ======
#----------------------------------------------------------------------#
emis<-NULL
for (x in 1:length(files)){
  emis[[x]]<-readEmissionTable(paste(modelDir, files[x], sep = '/'))%>%
    as.data.frame()
  tmp<-row.match(matrix(emis[[x]]), uniq)
  emis[[x]]<-as.matrix(emis[[x]])
  rownames(emis[[x]])<-tmp
}

emis

plots<-NULL
for (z in 1){
  euc<-matrix(data =NA, nrow=nrow(emis[[z]]), ncol=nrow(emis[[z+1]]))
  for (x in 1:nrow(emis[[z]])){
    for (y in 1:nrow(emis[[z+1]])){
      euc[x, y]<-euclidean(emis[[z]][x,], emis[[z+1]][y,])
    }
  }
  euc<-euc/max(euc)
  
  rownames(euc)<-rownames(emis[[z]])
  colnames(euc)<-rownames(emis[[z+1]])
  
  euc<-euc[mixedsort(rownames(euc)), mixedsort(colnames(euc))]
  
  plots[[z]]<-grid.grabExpr(draw(Heatmap(euc, col = brewer.pal(10,"YlOrRd"), show_row_dend = FALSE, 
                                         show_column_dend =FALSE,show_heatmap_legend = FALSE, 
                                         cluster_rows = FALSE, cluster_columns = FALSE)))
}
grid.newpage()
grid.arrange(grobs=plots, ncol=5, layout)


plots[[3]]<-grid.grabExpr(draw(Heatmap(euc, col = brewer.pal(10,"YlOrRd"), show_row_dend = FALSE, 
        show_column_dend =FALSE,show_heatmap_legend = FALSE, 
        cluster_rows = FALSE, cluster_columns = FALSE)))


df

#----------------------------------------------------------------------#
#             ==========================================================
#----------------------------------------------------------------------#
##try upsetr plot
## with rows as combinations of marks
## and showing how many times each is found

#----------------------------------------------------------------------#
# Genomic coverage of each   ================================
#----------------------------------------------------------------------#
## add the States in the order which they are found according to the
## master unique states
set.seed(150)
brbg <- c(brewer.pal(12, 'Set3'), brewer.pal(3, 'Dark2')) %>%
  sample()

cov<-matrix(NA, nrow=110, ncol = 3) #110, 12
cov[,1]<-unlist(lapply(emis, rownames))

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
cov

cov<-as.data.frame(cov) %>%
  setnames(c('State', 'Model', 'genomeCov'))
cov$genomeCov<-as.numeric(levels(cov$genomeCov))[cov$genomeCov]

ggplot(cov, aes(fill = State, x = Model, y = genomeCov))+
  geom_bar(stat = 'identity', position = 'stack', width = 0.7, colour = "black") +
  scale_fill_manual(values = brbg)+
  theme_bw()

#----------------------------------------------------------------------#
#             ==========================================================
#----------------------------------------------------------------------#
