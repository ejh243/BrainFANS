##---------------------------------------------------------------------#
##
## Title: Collate Summary Statistics ===================================
##
## Purpose of script: to collate the summary statistics of 
##                    the stage 1 QC metrics
##
## Author: Jessica Shields
##
## Date Created: 2022-07-06
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 
## set working directory
setwd("")

#----------------------------------------------------------------------#
# LOAD PACKAGES ========================================================
#----------------------------------------------------------------------#
project='epiGaba'
source("ChIPSeq/config/config.r")




library(ChIPQC)
library(knitr)
library(dplyr)
library(reshape)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel")
library(org.Hs.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
register(DoparParam())
registered() 
bpparam("SerialParam")

#----------------------------------------------------------------------#
# IMPORT DATA ==========================================================
#----------------------------------------------------------------------#

#setwd(dataDir)
sampleSheet<-read.csv(paste(metaDir, "sampleSheetForChipQC.csv",sep = "/"), stringsAsFactors = FALSE)

#if combined chipseqQC object does not exist, create it
  if (file.exists(paste0(peakDir, "/QCOutput/ChIPQCObject.rdata"))==FALSE){
  chipFiles <- list.files(paste0(peakDir, "/QCOutput"), pattern = "ChIPQCObject_")
  
  datALL<-NULL
  for (each in chipFiles) {
    print(each)
    load(paste(peakDir, "QCOutput", each, sep = "/"))
    datALL<- c(datALL, QCsample(dat))
  }
  # create combined chipqcexperiment object
  dat<- ChIPQC(sampleSheet, samples=datALL)
  save(dat, file = paste0(peakDir, "/QCOutput/ChIPQCObject.rdata"))
} else {
  load(paste0(peakDir, "/QCOutput/ChIPQCObject.rdata"))
}



##doesnt work
sumTab<-cbind(sampleSheet[,c(1,2)], reads(dat), unlist(lapply(peaks(dat), length)))
colnames(sumTab)<-c(colnames(sampleSheet)[c(1,2)], "Reads", "Peaks")
kable(sumTab, caption = "Table 1: Summary of aligned reads and number of peaks")



qcOut<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = ".pbc.qc")
sampleIDs<-gsub(".filt.srt.nodup.pbc.qc", "", qcOut)
qcDat<-NULL
for(each in qcOut){
  print(each)
  qcDat<-rbind(qcDat, read.table(paste0(alignedDir, "/ENCODEMetrics", "/", each)))
}
rownames(qcDat)<-sampleIDs
colnames(qcDat)<-c("TotalReadPairs","DistinctReadPairs","OneReadPair","TwoReadPairs","NRF","PBC1","PBC2")

## add ENCODE classification
qcDat$Complexity<-cut(qcDat$NRF, c(0,0.5,0.8,0.9, 1), labels = c("Concerning", "Acceptable", "Compliant", "Ideal"))
qcDat$BottleneckingLevel1<-cut(qcDat$PBC1, c(0,0.5,0.8,0.9, 1), labels = c("Severe", "Moderate", "Mild", "None"))
qcDat$BottleneckingLevel2<-cut(qcDat$PBC2, c(0, 1,3,10,100), labels = c("Severe", "Moderate", "Mild", "None"))

kable(qcDat[,c("Complexity", "BottleneckingLevel1", "BottleneckingLevel2")], caption = "Table 2: Summary of ENCODE quality control metrics")


## Enrichment of ChIP
par(mfrow = c(1,1))
hist(ssd(dat), cex.axis = 2, cex.lab = 2, xlab = "SSD", ylab = "nSamples", main = "")



## Specificity of reads in peaks
plotRap(dat, facet = FALSE) + theme(text = element_text(size=20), 
                                    axis.title = element_text(size = 10, angle = 25),
                                    axis.text.x = element_text(angle=90))


par(mar = c(5,5,1,1))
y_lim<-c(0,max(frip(dat)*100, na.rm = TRUE)*1.1)
barplot(frip(dat)*100, ylim = y_lim, ylab = "% Reads in Peaks", names.arg = sampleSheet$IID, 
        cex.names = 2, cex.axis = 2, cex.lab = 2)



## Relative enrichment in genomic intervals
plotRegi(dat) + theme(text = element_text(size=20))


## only works for samples with > 0 peaks
nPeaks<-unlist(lapply(peaks(dat), length))
peakAnnoList <- lapply(peaks(dat)[which(nPeaks > 0)], annotatePeak,tssRegion=c(-500, 500), 
                       TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                       annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnnoList) + theme(text = element_text(size=20))

plotDistToTSS(peakAnnoList) + theme(text = element_text(size=20))

## Sample clustering

plotCorHeatmap(dat) + theme(text = element_text(size=20))

plotPrincomp(dat, dotSize=2, cex.axis = 2, cex.lab = 2) 



## Pathway analysis of annotated genes

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

#----------------------------------------------------------------------#
# ALIGNMENT AND FILTERING STATS ========================================
#----------------------------------------------------------------------#
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)

## PROGRESS SUMMARY
processSum <- read.csv(paste0(metaDir, "/summariseSampleProcessingProgress.csv"), stringsAsFactors = FALSE, strip.white = TRUE)

## MULTIQC
fastqc<-read.table(paste0(fastQCDir, "/multiqc/multiqc_data/multiqc_fastqc.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ALIGNMENT STATISTICS
alignQC<-read.table(paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ENCODE METRICS
flagstat<-read.csv(paste0(alignedDir, "/ENCODEMetrics/collateFlagStatMetrics.txt"), header = FALSE)

## FRIP
fripFiles <- list.files(paste0(peakDir, "/QCOutput"), pattern = "FRIP")
#----------------------------------------------------------------------#
# WRANGLE DATA =========================================================
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

table(pheno$fraction)

## b. PROGRESS SUMMARY
processSum <- processSum[match(c(pheno$sampleID, unique(pheno$controlID)), processSum$sampleID),]

## c. MULTIQC
## wrangle multiqc stats to retain important columns
colsKeep<-c("Total.Sequences","total_deduplicated_percentage", "Sequences.flagged.as.poor.quality")
mergeStats<-cbind(fastqc[match(processSum$R1Filename, fastqc$Filename),colsKeep], fastqc[match(processSum$R2Filename, fastqc$Filename),colsKeep])

## d. MITOCHONDRIAL READS
#mergeStats<-cbind(mergeStats, readCounts$V2/readCounts$V3*100)
#colnames(mergeStats)[ncol(mergeStats)]<-"PercentMTReads"

## e. ALIGNMENT QC
aIndex<-match(processSum$sampleID, gsub("\\.bowtie", "", alignQC$Sample))
mergeStats<-cbind(mergeStats,alignQC[,c("overall_alignment_rate", "paired_total")])

## f. ENCODE METRICS
## load flagstat metrics calculated as part of encode qc pipeline
flagstat<-flagstat[match(processSum$sampleID, gsub('.filt.srt.nodup', '', flagstat$V1)),]
files<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = ".pbc.qc")
eMetrics<-NULL
for(each in files){
  tmp<-read.table(paste0(alignedDir, "/ENCODEMetrics/", each))
  eMetrics<-rbind(eMetrics, tmp)
}
eMetrics
colnames(eMetrics)<-c("TotalReadPairs","DistinctReadPairs","OneReadPair","TwoReadPairs","NRF","PBC1","PBC2")
eIndex<-match(processSum$sampleID, gsub(".pbc.qc", "", files))
eMetrics<-eMetrics[eIndex,]
mergeStats<-cbind(mergeStats,eMetrics)

## 6. FRIP
if ( length(fripFiles) != 0 ) {
  fripStats<-read.csv(paste0(peakDir, "/QCOutput/", fripFiles[1]))
  for(each in fripFiles[-1]){
    fripStats<-rbind(fripStats, read.csv(paste0(peakDir, "/QCOutput/", each)))
  }
  
  fripStats<-fripStats[match(processSum$sampleID, fripStats$SampleName),]
} else { warning('FRIP stats do not appear to have been calculated') }

fripStats$FripMACS2PE <- fripStats$ReadsinMACS2PEPeaks/fripStats$BAMTotalReads
#----------------------------------------------------------------------#
# COUNT SAMPLES ========================================================
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

## count number of samples with X million useable fragments
# 45 million usable fragments - broad
# 20 million usable fragments - narrow

broad<-c('H3F3A', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K79me2', 'H3K79me3', 'H3K9me1', 'H3K9me2', 'H4K20me1')
narrow<-c('H2AFZ', 'H3ac', 'H3K27ac', 'H3K4me2', 'H3K4me3', 'H3K9ac')

readThres<-seq(0,max(mergeStats[,1], na.rm = TRUE)+10^6, by = 10^6)
nSamples<-matrix(data = NA, nrow = length(readThres), ncol = 3)
colnames(nSamples)<-c('all', 'narrow', 'broad')

for(i in readThres){
  nSamples[1+(i/10^6),1] <- sum(mergeStats[,1] > i, na.rm = TRUE)
  colNum<-2
  nSamples[1+(i/10^6),colNum] <- (sum(mergeStats[which(c(pheno$target) %in% narrow),1] > i, na.rm = TRUE)*2)
  colNum<-colNum+1
  nSamples[1+(i/10^6),colNum] <- (sum(mergeStats[which(c(pheno$target) %in% broad),1] > i, na.rm = TRUE)*2)
  colNum<-colNum+1
}

#----------------------------------------------------------------------#
#
#----------------------------------------------------------------------#
par(mfrow = c(2,2))
hist(mergeStats[,1], main = "R1", xlab = "Number of Reads", ylab = "Number of Samples", breaks = 100)
hist(mergeStats[,4], main = "R2", xlab = "Number of Reads", ylab = "Number of Samples", breaks = 100)
hist(mergeStats[,2], main = "R1", xlab = "% unique reads", ylab = "Number of Samples", breaks = 100)
hist(mergeStats[,5], main = "R2", xlab = "% unique reads", ylab = "Number of Samples", breaks = 100)



## a. Plot alignment rate against number of reads
plot(mergeStats[,"paired_total"], mergeStats[,"overall_alignment_rate"], xlab = "Number of Reads", ylab = "Alignment rate", pch = 16, col = colorBlindGrey8[as.factor(pheno$fraction)])
abline(h = 95, lty = 2)
abline(h = 80, lty = 2, col = "#66C2A5")

## plot minimum number of useable fragments
#broad
par(mfrow = c(1,1))
plot(readThres/10^6, nSamples[,1], type = "l", lwd = 2, xlab = "Minimum number of reads (millions)", ylab = "Number of Samples")
axis(4, seq(0,100, 20), at = nrow(pheno)*seq(0,1,0.2))
mtext(side = 4, "% of samples", line = 2)
abline(v = (25)/0.8, lty = 2)
abline(v = (50), lty = 2)
abline(v = median(mergeStats[,1], na.rm = TRUE)/10^6, col = colorBlindGrey8[2], lty = 2)
# if there is more than one cohort
if ( ncol(nSamples) > 2 ) {
  plot(readThres/10^6, nSamples[,2], type = "l", lwd = 2, xlab = "Minimum number of reads (millions)", ylab = "Number of Samples")
  abline(v = (25)/0.8, lty = 2)
  abline(v = (50), lty = 2)
  abline(v = median(mergeStats[,1], na.rm = TRUE)/10^6, col = colorBlindGrey8[2], lty = 2)
  for(i in 3:ncol(nSamples)){
    lines(readThres/10^6, nSamples[,i], lwd = 2)
  }

  plot(readThres/10^6, nSamples[,2]/nSamples[1,2]*100, type = "l", lwd = 2, xlab = "Minimum number of reads (millions)", ylab = "% of Samples", ylim = c(0,100))
  abline(v = (25)/0.8, lty = 2)
  abline(v = (50), lty = 2)
  abline(v = median(mergeStats[,1], na.rm = TRUE)/10^6, col = colorBlindGrey8[2], lty = 2)
  for(i in 3:ncol(nSamples)){
    lines(readThres/10^6, nSamples[,i]/nSamples[1,i]*100, lwd = 2)
  }
}

## d. Plot %unique reads against read number
plot(mergeStats[,1], mergeStats[,2], xlab = "Number of Reads", ylab = "% unique reads", pch = 16, col = colorBlindGrey8[as.factor(pheno$fraction)])
legend('topright', legend = unique(pheno$fraction), 
       fill = unique(colorBlindGrey8[as.factor(pheno$fraction)]),cex=0.8) # does this always work??

#
plot(mergeStats[,1], mergeStats$DistinctReadPairs, xlab = "Number of reads", ylab = "Number of distinct aligned reads", pch = 16, col = colorBlindGrey8[as.factor(pheno$fraction)])
legend('topright', legend = unique(pheno$fraction), 
       fill = unique(colorBlindGrey8[as.factor(pheno$fraction)]),cex=0.8) # does this always work??

par(mfrow = c(3,1))
hist(eMetrics[,"NRF"], breaks = seq(min(c(eMetrics[,"NRF"], 0.65), na.rm = TRUE), max(c(eMetrics[,"NRF"], 0.95), na.rm = TRUE), length.out = 10), main = "", xlab = "Non-Redundant Fraction (NRF)", ylab = "nSamples")
abline(v = 0.7, lty = 2)
abline(v = 0.9, lty = 2)
mtext("Concerning", at = 0.65, side = 3)
mtext("Acceptable", at = 0.8, side = 3)
mtext("Ideal", at = 0.95, side = 3)
hist(eMetrics[,"PBC1"], breaks = seq(min(c(eMetrics[,"PBC1"], 0.65), na.rm = TRUE), max(c(eMetrics[,"PBC1"], 0.95), na.rm = TRUE), length.out = 10), main = "", xlab = "PCR Bottlenecking Coefficient 1 (PBC1)", ylab = "nSamples")
abline(v = 0.7, lty = 2)
abline(v = 0.9, lty = 2)
mtext("Severe", at = 0.65, side = 3)
mtext("Moderate", at = 0.8, side = 3)
mtext("None", at = 0.95, side = 3)
hist(eMetrics[,"PBC2"], breaks = seq(min(c(eMetrics[,"PBC2"], 0.95), na.rm = TRUE), max(c(eMetrics[,"PBC2"], 3.05), na.rm = TRUE), length.out = 10), main = "", xlab = "PCR Bottlenecking Coefficient 2 (PBC2)", ylab = "nSamples")
abline(v = 1, lty = 2)
abline(v = 3, lty = 2)
mtext("Severe", at = 0.95, side = 3)
mtext("Moderate", at = 2, side = 3)
mtext("None", at = 3.05, side = 3)

#----------------------------------------------------------------------#
# CHIPQC PLOTS ========================================================
#----------------------------------------------------------------------#

## Coverage histogram
plotCoverageHist(dat,colourBy=('Factor'), facetBy=("Tissue")) +
  theme_minimal() +
  scale_colour_brewer(palette = "Set2")

#plot ssd
ssdplot<-plotSSD(dat, facetBy=c("Tissue", "Factor"))+
  theme_minimal() +
  scale_colour_brewer(palette = "Set2")

ssdplot + scale_size(range=4:5)

## Cross-coverage
ccplot<-plotCC(dat, colourBy=('Factor'), facetBy=("Tissue"))+
  theme_minimal()+
  scale_y_continuous(limits=c(0,0.05))
ccplot+ scale_colour_brewer(palette = "Set2")
ccplot+facet_wrap(~Tissue,scales="free_y")


plotCC(dat, facetBy=c("Tissue","Condition"),colourBy = "Factor")+
  theme_minimal()


## Plotting Relative Enrichment of reads in Genomic Intervals
plotRegi(dat, facetBy=c("Tissue", "Factor")) + 
  theme_minimal()+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=90))+
  scale_fill_viridis_c()

# feature distribution
nPeaks<-unlist(lapply(peaks(dat), length))
peakAnnoList <- lapply(peaks(dat)[which(nPeaks > 0)], annotatePeak,tssRegion=c(-500, 500), 
                       TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                       annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnnoList) + theme(text = element_text(size=10))+
  theme_minimal()



##Plotting Peak Profiles
plotPeakProfile(dat, colourBy=('Factor'), facetBy=("Tissue"))+
  theme_minimal()+
  scale_colour_brewer(palette = "Set2")

# Plotting Reads overlapping Peaks and the Blacklist
## Specificity of reads in peaks
plotRap(dat, facet = FALSE) + 
  theme_minimal()+
  theme(text = element_text(size=10), 
        axis.title.y = element_text(size = 10, angle = 180),
        axis.text.x = element_text(angle=90)) +
  scale_fill_brewer(palette = "Set2")

plotFrip(dat, facet = FALSE) + 
  theme_minimal()+
  theme(text = element_text(size=10), 
        axis.title = element_text(size = 10, angle = 25),
        axis.text.x = element_text(angle=90))

plotFribl(dat, facet = FALSE) + 
  theme_minimal()+
  theme(text = element_text(size=10), 
        axis.title = element_text(size = 10, angle = 25),
        axis.text.x = element_text(angle=90))


##  Plotting Sample Clustering
plotCorHeatmap(dat,attributes=c("Tissue","Factor"))
plotPrincomp(dat,attributes=c("Tissue","Factor"))


