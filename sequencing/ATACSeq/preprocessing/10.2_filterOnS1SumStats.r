##---------------------------------------------------------------------#
##
## Title: Stage 1 Filter
##
## Purpose of script: to reformat multiQC output into more meaningful summaries
##
## Author: Eilis Hannon, edited by Jessica Shields
##
## Date Created: 2022-02-11
##
##---------------------------------------------------------------------#
## clear the R environment
rm(list=ls())

## load arguments
args = commandArgs(trailingOnly=TRUE)

project<-args[1]


## set default if no args supplied
args= as.numeric(args)
if (length(args)==1) {
  warning("No filtering parameters specified, using default parameters:
          Number of reads > 10 million
          Alignment rate > 80%
          Number of filtered/aligned reads > 20 million")
  args[2] = 10 #read number
  args[3] = 80 # alignment rate
  args[4] = 20 # aligned/filtered
}

source("ATACSeq/config/config.r")

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(vioplot)
library(FME)


## create colourblind friendly palette
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
# METADATA
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)

## PROGRESS SUMMARY
processSum <- read.csv(paste0(metaDir, "/summariseSampleProcessingProgress.csv"), stringsAsFactors = FALSE, strip.white = TRUE)

## MULTIQC
fastqc<-read.table(paste0(fastQCDir, "/multiqc/multiqc_data/multiqc_fastqc.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## MITOCHONDRIAL READS
readCounts<-read.table(paste0(alignedDir, "/countMTReads.txt"), fill = TRUE, skip = 1)

## ALIGNMENT STATISTICS
alignQC<-read.table(paste0(alignedDir, "/multiqc/multiqc_data/multiqc_bowtie2.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ENCODE METRICS
flagstat<-read.csv(paste0(alignedDir, "/ENCODEMetrics/collateFlagStatMetrics.txt"), header = FALSE)

## FRIP
fripFiles <- list.files(paste0(peakDir, "/QCOutput"), pattern = "FRIP")

#----------------------------------------------------------------------#
# WRANGLE DATA
#----------------------------------------------------------------------#
## a. METADATA
## check metadata is in the correct format

# METADATA REQUIREMENTS:
# essential columns: "sampleID",	"cohort",	"fraction",	"experiment",	"individualID"
# samples can be identified as unique from sampleID, cohort, fraction and individualID

missingCol<-setdiff(c("sequencingBatch", "sampleID",	"cohort",	"fraction",	"experiment",	"individualID"), colnames(pheno))
if (length(missingCol) != 0 ) {
  if ("sequencingBatch" %in% missingCol){
    pheno$sequencingBatch<-pheno$cohort
  } else {
    warning('Incorrect column names. Missing column(s): ', missingCol)
  }
}

## exclude duplicates from metadata
dups<-names(which(table(paste(pheno$individualID, pheno$fraction, pheno$tissue, sep = "_")) > 1))
pheno<-unique(pheno[,intersect(c("sequencingBatch", "sampleID","cohort","fraction","individualID", "tissue"), colnames(pheno))])
pheno$sequencingRuns<-1
pheno$sequencingRuns[match(dups, paste(pheno$individualID, pheno$fraction, pheno$tissue, sep = "_"))]<-2
#sequencingRuns is not used anywhere else?? 

table(pheno$fraction)
barplot(table(table(pheno$individualID)), xlab = "Number of Samples", ylab = 'Number of individuals', col="white")

## b. PROGRESS SUMMARY
## sort if processSum$sampleID is not equivalent to pheno$sampleID

processSum$sampleName<-gsub("^[[:digit:]]...._", "", gsub("_S[[:digit:]]+$", "", processSum$sampleID))
index<-match(pheno$sampleID, processSum$sampleName)
index[is.na(index)]<-match(pheno$sampleID[is.na(index)], processSum$sampleName)
processSum <- processSum[index,]

## c. MULTIQC
## wrangle multiqc stats to retain important columns

colsKeep<-c("Total.Sequences","total_deduplicated_percentage", "Sequences.flagged.as.poor.quality")
mergeStats<-cbind(fastqc[match(processSum$R1Filename, fastqc$Filename),colsKeep], fastqc[match(processSum$R2Filename, fastqc$Filename),colsKeep])

## d. MITOCHONDRIAL READS

readCountsMT<-readCounts[grep("_postFilter_statsperchr.txt", readCounts$V1, invert = TRUE),]
readCountsMT$V1<-gsub("_statsperchr.txt", "", readCountsMT$V1)

readCountsMT<-readCountsMT[match(processSum$sampleID, readCountsMT$V1),]
mergeStats<-cbind(mergeStats, readCountsMT$V2/readCountsMT$V3*100)
colnames(mergeStats)[ncol(mergeStats)]<-"PercentMTReads"

## e. ALIGNMENT QC

aIndex<-match(processSum$sampleName, gsub("^[[:digit:]]...._", "", gsub("\\.bowtie", "", alignQC$Sample)))
mergeStats<-cbind(mergeStats,alignQC[aIndex,c("overall_alignment_rate", "paired_total")])

## f. ENCODE METRICS
## load flagstat metrics calculated as part of encode qc pipeline

flagstat<-flagstat[match(processSum$sampleID, flagstat$V1),]

files<-list.files(paste0(alignedDir, "/ENCODEMetrics"), pattern = ".pbc.qc")
eMetrics<-NULL
for(each in files){
  tmp<-read.table(paste0(alignedDir, "/ENCODEMetrics/", each))
  eMetrics<-rbind(eMetrics, tmp)
}
eMetrics
colnames(eMetrics)<-c("TotalReadPairs","DistinctReadPairs","OneReadPair","TwoReadPairs","NRF","PBC1","PBC2")

#eIndex<-match(processSum$sampleID, gsub("_sorted_chr1.bam.nodup.pbc.qc", "", files))
#eMetrics<-eMetrics[eIndex,]

mergeStats<-cbind(mergeStats,eMetrics)

## 6. FRIP
if ( length(fripFiles) != 0 ) {
  fripStats<-read.csv(paste0(peakDir, "/QCOutput/", fripFiles[1]))
  for(each in fripFiles[-1]){
    fripStats<-rbind(fripStats, read.csv(paste0(peakDir, "/QCOutput/", each)))
  }
  
  fripStats<-fripStats[match(processSum$sampleID, fripStats$SampleName),]
} else { warning('FRIP stats do not appear to have been calculated') }

fripStats$FripMACS2TagAlign <- fripStats$MACS2TagAlignPeaks/fripStats$BAMTotalReads
fripStats$FripMACS2PE <- fripStats$MACS2PEPeaks/fripStats$BAMTotalReads

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
# PLOT NUMBER OF READS
#----------------------------------------------------------------------#
## a. Plot distribution of number of sequencing reads

## create summary table
statTable<-table(pheno$sequencingBatch)
statTable<-cbind(statTable, aggregate(mergeStats[,1], by = list(pheno$sequencingBatch), quantile, c(0.25,0.5,0.75), na.rm = TRUE))
statTable<-cbind(statTable, aggregate(mergeStats[,2], by = list(pheno$sequencingBatch), quantile, c(0.25,0.5,0.75), na.rm = TRUE)[,-1])
statTable<-cbind(statTable, aggregate(mergeStats[,1] < 10*10^6, by = list(pheno$sequencingBatch), sum, na.rm = TRUE)[,-1])
statTable$PercentageLessThan10Million<-statTable[,ncol(statTable)]/statTable[,2]*100
write.csv(statTable, file = paste0(qcDir, "/summarySequencingReadsByProject.csv"))

## create object for samples not passing the min read threshold
lowReadCounts<-pheno[which(mergeStats[,1] < args[2]*10^6),]

## create object for samples not passing the min alignment threshold

lowAlignmentRate<-pheno[which(mergeStats[,"overall_alignment_rate"] < args[3]),]

#----------------------------------------------------------------------#
# PERIODICITY
#----------------------------------------------------------------------#
## load ATACQC metrics 
propNucleosomesAll<-NULL
filePaths <- list.files(paste0(alignedDir,"/QCOutput/"), pattern = "FragmentDistribution_Batch")
for(each in filePaths){
  load(paste0(alignedDir,"/QCOutput/", each))
  propNucleosomesAll<-rbind(propNucleosomesAll, propNucleosomes)
}

propNucleosomesAll<-propNucleosomesAll[match(processSum$sampleID, rownames(propNucleosomesAll)),]

## a successful ATAC experiment should have periodicity in the fragment distribution
## detect samples with no evidence of mono, di, tri nucleosomes
propNucleosomesAll > 0.9

## a successful ATAC experiment should have decreasing proportions in the nucleosome, mono, di,tri- nucleosomes
## define a tolerence for these comparisions
tol <- 0.05
decreasingProps<-cbind(propNucleosomesAll[,1]+tol > propNucleosomesAll[,2], 
                       propNucleosomesAll[,2]+tol > propNucleosomesAll[,3], 
                       propNucleosomesAll[,3]+tol > propNucleosomesAll[,4]
)
colnames(decreasingProps)<-c("NFR>Mono", "Mono>Di", "Di>Tri")

#----------------------------------------------------------------------#
# COLLATE FRIP AND PEAK STATISTICS
#----------------------------------------------------------------------#
#if ()
### summarise number of aligned, deduplicated reads
readThres<-seq(0,max(fripStats$BAMTotalReads, na.rm = TRUE)+10^6, 10^6)
nSamples<-matrix(data = NA, nrow = length(readThres), ncol = length(unique(pheno$sequencingBatch))+1)
for(i in readThres){
  nSamples[1+(i/10^6),1] <- sum(fripStats$BAMTotalReads > i, na.rm = TRUE)
  colNum<-2
  for(each in unique(pheno$sequencingBatch)){
    nSamples[1+(i/10^6),colNum] <- sum(fripStats$BAMTotalReads[which(pheno$sequencingBatch == each)] > i, na.rm = TRUE)
    colNum<-colNum+1
  }
}

## calculate frip
fripStats$FripMACS2TagAlign <- fripStats$ReadsinMACS2TagAlignPeaks/fripStats$BAMTotalReads
fripStats$FripMACS2PE <- fripStats$ReadsinMACS2PEPeaks/fripStats$BAMTotalReads

#----------------------------------------------------------------------#
# FILTER SAMPLES
#----------------------------------------------------------------------#

QCPASS<-cbind(mergeStats$overall_alignment_rate > args[3], #alignment rate 
              fripStats$BAMTotalReads > args[4]*10^6, #total filtered/aligned reads
              eMetrics$NRF > 0.7, #encode
              eMetrics$PBC1 > 0.7, #encode
              eMetrics$PBC2 > 1, #encode
              propNucleosomesAll[,2] > 0.15, #>15% fractions mononucleosome 
              decreasingProps[,1]) #periodicity

keep<-rowSums(QCPASS) == ncol(QCPASS)
length(keep[keep==TRUE])

write.table(names(keep[keep==TRUE]), file = paste0(metaDir, "/stage1QCSamples.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)




