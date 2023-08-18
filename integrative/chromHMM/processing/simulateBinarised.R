##---------------------------------------------------------------------#
##
## Title: Simulate chromHMM data =======================================
##
## Purpose of script: to simulate binary matrix for chromHMM input 
##
## Author: Jessica Shields
##
## Date Created: 2022-10-17
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 
## set working directory
setwd("/lustre/home/jms260/BrainFANS")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args<-"roadMapPFC"
} 
intproject<-args
simproject<-'simulateRM'
type <- 'PFC' #'N+'

source("integrative/chromHMM/config/config.r")

# create pathway for simulated data
simulateDir<-gsub("[^/]+$", paste0(simproject, "/2_mergeBinarised"), chromDir) # catch everything after the last "/"

dir.create(simulateDir, recursive=TRUE)

#----------------------------------------------------------------------#
# LOAD PACKAGES AND FUNCTIONS ==========================================
#----------------------------------------------------------------------#
library(GenomicFeatures)
library(magrittr)

#----------------------------------------------------------------------#
# Generate chromosome length vector and proportions of combinations ====
#----------------------------------------------------------------------#
# vector for chr1:22
chr<-paste0('chr', seq(1:22)) %>% sort()

# vector for chromosome lengths (divided by chromhmm 200bp bins)
chrLength<-getChromInfoFromUCSC("hg38") %>%
  .[1:24,2] %>%
  sapply(FUN = function(x) {round(x/200)})

# give chromosome number to length and order as character
chrOrder<-cbind(as.character(1:22), chrLength[1:22]) %>%
  .[order(.[,1]),] 

# table of 10 unique mark combinations
freqTable<-read.table(paste0(simulateDir, "/freq.txt"), header=TRUE)

freqTable<-
  list.files(modelDir, pattern = 'emissions_15') %>% 
    paste(modelDir, ., sep='/') %>% 
    readEmissionTable()

print(paste(nrow(unique(freqTable)), " true unique states"))

#----------------------------------------------------------------------#
# SIMULATE DATA =======================================================
#----------------------------------------------------------------------#
windowFile <- list.files(qcDir, pattern = 'dense.200.bed')
sim <- read.table(paste(qcDir, windowFile[1], sep = '/'))

df<-NULL
for (i in chr){
#  tmp <- sim %>% filter(V1==i)
#  df <- cbind(df, table(tmp$V2))
  tmp <- filter(sim, sim$V1==i) %>% 
    {freqTable[match(.[,2], rownames(freqTable)),1:(ncol(freqTable))]} %>%
    na.omit() %>%
    rbind(c(type, i, rep('', ncol(.)-2)), colnames(.), .)
#  write.table(tmp, paste0(simulateDir,"/", type, "_", i,"_binary.txt"), sep='\t', col.names = FALSE ,row.names = FALSE, quote = FALSE)
}


sim<-lapply(chr, function(i){
  filter(sim, V1==i) %>% 
    {freqTable[match(.[,2], rownames(freqTable)),1:(ncol(freqTable))]} %>%
    na.omit() %>%
    rbind(c(type, i, rep('', ncol(.)-2)), colnames(.), .) %>%
    write.table(paste0(simulateDir,"/", type, "_", i,"_binary.txt"), sep='\t', col.names = FALSE ,row.names = FALSE, quote = FALSE)
})
#----------------------------------------------------------------------#
# CHECKS ===============================================================
#----------------------------------------------------------------------#
# ensure correct states are present across the chromosomes
reshape2::melt(df) %>% ggplot(aes(Var2, value, fill=Var1))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()


# ensure that the created table is identical to the written.
tmp<-filter(sim, V1==chr[1]) %>% 
  {freqTable[match(.[,2], rownames(freqTable)),1:(ncol(freqTable)-1)]} %>%
  na.omit()
rownames(tmp)<-NULL
check1<-read.table(paste0(simulateDir, "/", type, "_", chr[1], "_binary.txt"), skip = 1, header = TRUE) %>% as.matrix()

identical(check1, tmp)

if (check1){
  print("Data is identical")
}

# ensure that there are no additional combinations
check2<-NULL
for (i in 1:22){
  tmp<-read.table(paste0(simulateDir, "/", type, "_chr", i, "_binary.txt"), skip = 1, header = TRUE) %>% unique()
  check2<-rbind(check2, tmp) %>% unique()
  print(i)
}

if (nrow(check2) == nrow(freqTable)){
  print('Correct number of unique combinations.')
} else {
  print(paste("Incorrect number of combinations: data has",nrow(check2)-nrow(freqTable), "difference"))
}

#----------------------------------------------------------------------#
# SIMULATE NOISE =======================================================
#----------------------------------------------------------------------#
## varying degrees of noise
pc<-c(2, 10, 20)
pc<-2
intproject<-paste0('simulate_', pc/2)

source("integrative/chromHMM/config/config.r")

dir.create(chromDir)
dir.create(mergeDir)

# set prgress bar
progressBar<-txtProgressBar()

set.seed(100)
simMarks1<-NULL
for (x in 1:22){ #numbered chromosomes
  set.seed(100)
  sim<-simulateState(freq, chrLength[x])
  set.seed(100)
  sim0<-simulateState(freq, chrLength[x])
  onePC<-round(nrow(sim)*(pc/100))
#  print(onePC)
  dif<-NULL
  for (y in 1:ncol(sim)){
    set.seed(100)
    sim[sample(nrow(sim), onePC), y]<-rbinom(onePC,1,(1-mean(sim[,y])))
    dif[[y]]<-length(which(sim0[,y] != sim[,y]))/length(sim0[,y])
  }
  simMarks[[x]]<-dif
#  break
#  sim<-rbind(c('N+', paste0('chr',x), '', ''), colnames(sim),sim)
#  write.table(sim, paste0(mergeDir,"/N+_chr",x,"_binary.txt"), sep='\t', col.names = FALSE ,row.names = FALSE, quote = FALSE)
  
  setTxtProgressBar(progressBar, x/22)
}
simMarks

# write the percetnage difference from no-noise data to text file
lapply(simMarks, mean) %>%
         unlist() %>%
         mean() %>%
  write.table(paste0(mergeDir, "/noise.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)



apply(sim[3:nrow(sim), ], 2, sum)
as.data.frame(sim) %>% group_by_all() %>% count

## write some checks for this
setdiff(sim0[,1],sim[,1])


a<-c(0,0,1,1,0)
b<-c(0,1,0,1,0)
length(which(a != b))/length(b)

#----------------------------------------------------------------------#
# PLOT DATA ============================================================
#----------------------------------------------------------------------#

## aic

dirs<-c('simulate', 'simulate_1', 'simulate_5') #, 'simulate_10', 'simulate_50')

plots<-NULL
x=1
for (intproject in dirs){
  source("integrative/chromHMM/config/config.r")
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
  
  plots[[x]]<-ggplot(likeli, aes(states, aic))+
    geom_point()+
    theme_bw()+
    geom_line()
  x<-x+1
}

grid.arrange(grobs=plots, ncol=2)
