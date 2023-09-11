##---------------------------------------------------------------------#
##
## Title: Collate Stage 1 Summary Statistics
##
## Purpose of script:
##
## Author: Jessica Shields
##
## Date Created: 2022-05-26
##
##---------------------------------------------------------------------#
##
## set working directory
setwd("")
## clear the R environment
rm(list=ls()) 

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)
library(plyr)
library(dplyr)
library(corrplot)


project<-"WGBS/rizzardi"
source("/lustre/projects/Research_Project-MRC190311/scripts/sequencing/BSSeq/config/config.r")

## create colourblind friendly palette
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
sampleid<-read.table(paste0(metaDir, '/samples.txt'))[,1]

# METADATA
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
## PROGRESS SUMMARY
processSum <- read.csv(file.path(metaDir, "/summariseSampleProcessingProgress.csv"), stringsAsFactors = FALSE, strip.white = TRUE)

## MULTIQC
fastqc<-read.table(file.path(fastQCDir, "/multiqc/multiqc_data/multiqc_fastqc.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ALIGNMENT STATISTICS
alignQC<-read.table(file.path(alignedDir, "/multiqc/multiqc_data/multiqc_bismark_alignment.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ENCODE METRICS
encode<-read.csv(file.path(alignedDir, "/ENCODEMetrics/collateENCODEQCMetrics.txt"), header = TRUE)

#----------------------------------------------------------------------#
# WRANGLE DATA
#----------------------------------------------------------------------#
## a. METADATA
## check metadata is in the correct format
# METADATA REQUIREMENTS:
# essential columns: "sampleID",  "cohort", "fraction", "experiment", "individualID"
# samples can be identified as unique from sampleID, cohort, fraction and individualID
missingCol<-setdiff(c("sequencingBatch", "sampleID",  "cohort", "fraction", "experiment", "individualID"), colnames(pheno))
if (length(missingCol) != 0 ) {
  if ("sequencingBatch" %in% missingCol){
    pheno$sequencingBatch<-pheno$cohort
    if (length(missingCol) > 1) {
      warning('Incorrect column names. Missing column(s): ', list(missingCol))
      stop()
    }
  }
}

pheno<- pheno[match(pheno$sampleID, sampleid),]
## exclude duplicates from metadata
dups<-names(which(table(paste(pheno$individualID, pheno$fraction, pheno$tissue, sep = "_")) > 1))
pheno<-unique(pheno[,intersect(c("sequencingBatch", "sampleID","cohort","fraction","individualID", "tissue"), colnames(pheno))])

table(pheno$fraction)


## c. MULTIQC
## wrangle multiqc stats to retain important columns
colsKeep<-c("Total.Sequences","total_deduplicated_percentage", "Sequences.flagged.as.poor.quality")
mergeStats<-cbind(fastqc[match(processSum$R1Filename, fastqc$Filename),colsKeep], fastqc[match(processSum$R2Filename, fastqc$Filename),colsKeep])


## e. ALIGNMENT QC
aIndex<-match(processSum$sampleID, gsub("\\_1_val_1", "", alignQC$Sample))
mergeStats<-cbind(mergeStats,alignQC[aIndex,c("percent_aligned", "total_reads", "aligned_reads")])

## f. ENCODE METRICS
## load flagstat metrics calculated as part of encode qc pipeline
encode<-encode[match(processSum$sampleID, encode$sampleID),]
files<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = ".qc")
eMetrics<-NULL
for(each in files){
  tmp<-read.table(paste0(alignedDir, "/ENCODEMetrics/", each))
  eMetrics<-rbind(eMetrics, t(tmp))
}

colnames(eMetrics)<-c("coverage","conversionEfficiency")
eIndex<-match(processSum$sampleID, gsub(".qc", "", files))
eMetrics<-eMetrics[eIndex,]
mergeStats<-cbind(mergeStats,eMetrics)

## CORRELATION
files<-list.files(paste0(methylDir, '/QCOutput/'), pattern = '.corr.qc')
files <- files[gsub('_', ' ', gsub('.corr.qc', '', files)) %in% pheno$fraction]
#files <- files[!(gsub('_', ' ', gsub('.corr.qc', '', files)) %in% pheno$fraction)]

cMetrics<-NULL
for (x in 1:length(files)){
  tmp<-read.table(paste0(methylDir, '/QCOutput/',files[x]), sep='\t') %>%
    as.matrix()
  V<-c(tmp)%>% 
    t() %>% 
    as.data.frame()
  cMetrics[[x]]<-V
}
cMetrics<-do.call(rbind.fill, cMetrics) %>% 
  t() %>% 
  as.data.frame() 

colnames(cMetrics)<-gsub(".corr.qc", "", files)
cMetrics<-pivot_longer(cMetrics, cols=1:ncol(cMetrics), names_to = 'tissue', values_to = 'corr') %>%
  drop_na()

#----------------------------------------------------------------------#
# COUNT SAMPLES
#----------------------------------------------------------------------#
## count number of samples with X million reads
readThres<-seq(0,max(mergeStats[,1], na.rm = TRUE)+10^6, 10^6)
nSamples<-matrix(data = NA, nrow = length(readThres), ncol = length(unique(pheno$sequencingBatch))+1)
for(i in readThres){
  nSamples[1+(i/10^6),1] <- sum(mergeStats[,1] > i, na.rm = TRUE)
  colNum<-2
  for(each in unique(pheno$sequencingBatch)){
    nSamples[1+(i/10^6),colNum] <- sum(mergeStats[which(pheno$sequencingBatch == each),1] > i, na.rm = TRUE)
    colNum<-colNum+1
  }
}

#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#
#To provide an overall summary here is a table of the number of raw reads, aligned reads (post filtering) and peaks called for each sample.



# plot histograms
xlab <- c('Number of Reads', '% unique reads' ,'' , 'Number of Reads', '% unique reads')
read<-c('R1','R1','','R2', 'R2')
p=NULL
for (x in 1:5){
  p[[x]]<-mergeStats %>% dplyr::select(x)%>%
          ggplot(., aes_string(x = names(.)[1]))+ # select x values based on column index rather than name
          geom_histogram(bins=100, fill ='white', colour = 'black')+
          labs(title = read[x], x = xlab[x], y = 'Number of Samples')+
          theme_bw()
}
p<-p[- 3] # remove the null list item
grid.arrange(grobs = p,ncol=2)


# plot minimum number of reads
coeff=nSamples[1,1]/100
data.frame(readThres/10^6, nSamples) %>% 
  melt(., id='readThres.10.6') %>%
  ggplot(., aes_string(x = names(.)[1], y = names(.)[3], colour = names(.)[2]))+
  geom_line() +
  labs(x = 'Minimum number of reads (millions)', y = 'Number of Samples')+
  theme_bw()+
  geom_vline(xintercept=c((25)/0.8, 50, median(mergeStats[,1], na.rm = TRUE)/10^6), 
             colour = c('black', 'black', colorBlindGrey8[2]),
             alpha=0.9, linetype="dashed", show.legend =FALSE)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~./coeff, name = '% of samples'))+
  theme(legend.position = "none")+
  scale_color_manual(values = rep('black', ncol(nSamples)))

## a. Plot alignment rate against number of reads
ggplot(mergeStats[,c(7,8)], aes(x = total_reads, y = percent_aligned, colour = pheno$fraction))+
  geom_point() +
  geom_hline(yintercept=c(95, 80), 
             colour = c('black', colorBlindGrey8[2]),
             alpha=0.9, linetype="dashed", show.legend =FALSE)+
  labs(x = 'Number of Reads', y = 'Alignment rate')+
  theme_bw() +
  scale_color_manual(values = colorBlindGrey8[1:3])

## b. Plot alignment rate
p1<-data.frame(pheno$fraction, mergeStats$percent_aligned) %>%
  ggplot(., aes(x=pheno.fraction, y=mergeStats.percent_aligned, ))+
    theme_bw()+
    labs(x = 'Cell fraction', y = 'Alignment rate')+
    geom_violin(fill = "grey80")+
    geom_boxplot(width=0.05, fill='black')+
    stat_summary(fun.data=mean_sdl, mult=1, 
                 geom="point", color="white")

p2<-data.frame(pheno$fraction, mergeStats$aligned_reads) %>%
  ggplot(., aes(x=pheno.fraction, y=mergeStats.aligned_reads, ))+
  theme_bw()+
  labs(x = 'Cell fraction', y = 'Number of aligned reads')+
  geom_violin(fill = "grey80")+
  geom_boxplot(width=0.05, fill='black')+
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="point", color="white")

grid.arrange(grobs = list(p1, p2), nrow = 1)


plot(mergeStats[,1], mergeStats$aligned_reads, xlab = "Number of reads", ylab = "Number of distinct aligned reads", pch = 16, col = colorBlindGrey8[as.factor(pheno$fraction)])

data.frame(mergeStats[,1], mergeStats$aligned_reads, pheno$fraction) %>% 
    ggplot(., aes_string(x = names(.)[1], y = names(.)[2], colour = names(.)[3]))+
  geom_point() +
  labs(x = 'Number of reads', y = 'Number of distinct aligned reads')+
  theme_bw()+
  scale_colour_manual(values = colorBlindGrey8[1:3])+
  theme(legend.title = element_blank())

#coverage
data.frame(mergeStats[,1],mergeStats$coverage, pheno$fraction) %>%
  ggplot(., aes_string(x = names(.)[3], y = names(.)[2]))+
  theme_bw()+
  theme(axis.text = element_text(size = 10))+
  labs(x = 'Cell fraction', y = 'Coverage')+
  geom_violin(fill = "grey80")+
  geom_boxplot(width=0.05, fill='black')+
  stat_summary(fun.data=mean_sdl,  
               geom="point", color="white")+
  geom_hline(yintercept=30, alpha=0.7, linetype="dashed", show.legend =FALSE)

#----------------------------------------------------------------------#
# PLOT CONVERSION EFFICIENCY
#----------------------------------------------------------------------#
encode<-pivot_longer(encode, cols = 3:4)

data.frame(mergeStats$conversionEfficiency, pheno$sequencingBatch) %>%
  ggplot(., aes(x=pheno.sequencingBatch,y=mergeStats.conversionEfficiency))+
    theme_bw()+
    theme(axis.text = element_text(size = 10))+
    labs(x = '', y = 'Conversion efficiency %')+
    geom_violin(fill = "grey80")+
    geom_boxplot(width=0.05, fill='black')+
    stat_summary(fun.data=mean_sdl,
                 geom="point", color="white")+
    geom_hline(yintercept=98, alpha=0.7, linetype="dashed", show.legend =FALSE)


#----------------------------------------------------------------------#
# PLOT CORRELATION
#----------------------------------------------------------------------#
ggplot(cMetrics, aes(x=tissue, y=corr))+
  theme_bw()+
  theme(axis.text = element_text(size = 10))+
  xlab('')+
  ylim(c(0,1))+
  ylab('Pearson correlation')+
  geom_violin(fill = "grey80")+
  geom_boxplot(width=0.05, fill='black')+
  stat_summary(fun.data=mean_sdl, mult=1, 
               geom="point", color="white")+
  geom_hline(yintercept=0.8, alpha=0.7, linetype="dashed", show.legend =FALSE)


## create correlation plot 
corMergeStats<-cbind(rowMeans(mergeStats[,c(1,4)]),rowMeans(mergeStats[,c(2,5)]), rowMeans(mergeStats[,c(3,6)]),
                     mergeStats$percent_aligned,
                     mergeStats$aligned_reads, mergeStats[,c((ncol(mergeStats)-1):ncol(mergeStats))]) 
colnames(corMergeStats)<- c('total\nreads', 'dedup', 'poor\nqual', 'alignt\nrate', 'distinct\nreads', 'coverage', 'convrsn\nefficncy') 

corMergeStats<- corMergeStats[, colSums(corMergeStats, na.rm = TRUE) != 0]
corrplot.mixed(cor(corMergeStats, use = "p"), order = 'AOE', tl.cex=0.9, tl.col = 'black', upper.col = COL2('PiYG'), lower.col = COL2('PiYG'))


