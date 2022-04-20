##---------------------------------------------------------------------#
##
## Title: 
##
## Purpose of script:
##
## Author: 
##
## Date Created: 2022-03-25
##
##---------------------------------------------------------------------#

## clear the R environment
rm(list=ls()) 
## set working directory
setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/scripts/sequencing")


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
project='adChip'
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
# IMPORT DATA
#----------------------------------------------------------------------#

#setwd(dataDir)
sampleSheet<-read.csv(paste(metaDir, "sampleSheetForChipQC.csv",sep = "/"), stringsAsFactors = FALSE, row.names = 1)

#if combined chipseqQC object does not exist, create it
  if (file.exists(paste0(peakDir, "/QCOutput/ChIPQCObject.rdata"))==FALSE){
  chipFiles <- list.files(paste0(peakDir, "/QCOutput"), pattern = "ChIPQCObject_")
  
  datALL<-NULL
  for (each in chipFiles) {
    print(each)
    load(paste(peakDir, "QCOutput", each, sep = "/"))
    datALL<- c(datALL, QCsample(dat))
  }
  # Create combined chipqcexperiment object
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
# ALIGNMENT AND FILTERING STATS
#----------------------------------------------------------------------#

pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
## PROGRESS SUMMARY
processSum <- read.csv(paste0(metaDir, "/summariseSampleProcessingProgress.csv"), stringsAsFactors = FALSE, strip.white = TRUE)
## MULTIQC
fastqc<-read.table(paste0(fastQCDir, "/multiqc/multiqc_data/multiqc_fastqc.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
## MITOCHONDRIAL READS
#readCounts<-read.table(paste0(alignedDir, "/countMTReads.txt"), fill = TRUE, skip = 1)
## ALIGNMENT STATISTICS
alignQC<-read.table(paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
## ENCODE METRICS
#flagstat<-read.csv(paste0(alignedDir, "/ENCODEMetrics/collateFlagStatMetrics.txt"), header = FALSE)
## FRIP
fripFiles <- list.files(paste0(peakDir, "/QCOutput"), pattern = "FRIP")
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

table(pheno$fraction)

## b. PROGRESS SUMMARY


## c. MULTIQC
## wrangle multiqc stats to retain important columns
colsKeep<-c("Total.Sequences","total_deduplicated_percentage", "Sequences.flagged.as.poor.quality")
mergeStats<-cbind(fastqc[match(processSum$R1Filename, fastqc$Filename),colsKeep], fastqc[match(processSum$R2Filename, fastqc$Filename),colsKeep])
## d. MITOCHONDRIAL READS
#mergeStats<-cbind(mergeStats, readCounts$V2/readCounts$V3*100)
#colnames(mergeStats)[ncol(mergeStats)]<-"PercentMTReads"
## e. ALIGNMENT QC
#aIndex<-match(processSum$sampleName, gsub("^[[:digit:]]...._", "", gsub("_S[[:digit:]]+$", "", gsub("\\.bowtie", "", alignQC$Sample))))
mergeStats<-cbind(mergeStats,alignQC[,c("overall_alignment_rate", "paired_total")])
## f. ENCODE METRICS
## load flagstat metrics calculated as part of encode qc pipeline
#flagstat<-flagstat[match(processSum$sampleID, gsub('.filt.srt.nodup', '', flagstat$V1)),]
files<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = ".pbc.qc")
eMetrics<-NULL
for(each in files){
  tmp<-read.table(paste0(alignedDir, "/ENCODEMetrics/", each))
  eMetrics<-rbind(eMetrics, tmp)
}
eMetrics
colnames(eMetrics)<-c("TotalReadPairs","DistinctReadPairs","OneReadPair","TwoReadPairs","NRF","PBC1","PBC2")
eIndex<-match(processSum$sampleID, gsub(".filt.srt.nodup.pbc.qc", "", files))
eMetrics<-eMetrics[eIndex,]
mergeStats<-cbind(mergeStats,eMetrics)

## 6. FRIP
#if ( length(fripFiles) != 0 ) {
#  fripStats<-read.csv(paste0(peakDir, "/QCOutput/", fripFiles[1]))
#  for(each in fripFiles[-1]){
#    fripStats<-rbind(fripStats, read.csv(paste0(peakDir, "/QCOutput/", each)))
#  }
#  
#  fripStats<-fripStats[match(processSum$sampleID, fripStats$SampleName),]
#} else { warning('FRIP stats do not appear to have been calculated') }
#fripStats$FripMACS2TagAlign <- fripStats$MACS2TagAlignPeaks/fripStats$BAMTotalReads
#fripStats$FripMACS2PE <- fripStats$MACS2PEPeaks/fripStats$BAMTotalReads
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
# ChIPQC plots
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
  theme_minimal()
ccplot+ scale_colour_brewer(palette = "Set2")
ccplot+facet_wrap(~Tissue,scales="free_y")


plotCC(dat)+
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


