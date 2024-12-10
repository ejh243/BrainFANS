## ================================================================================================##
##                          ATAC-seq pipeline STEP 8.1: Nomalise counts                            ##
## ================================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/normCounts.r <project>                         ||
## - execute from scripts directory                                                                ||
##                                                                                                 ||
## DESCRIPTION: This scripts performs variance partition analysis (VPA) on the raw read counts     ||
## obtained in STEP 7.3 in peaks annotated as promoters. Then, it normalises these counts using    ||
## the RUV-III-NB method. VPA is performed again in the normalised counts.                         ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <project> : project on which analysis is being run                                            ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## -  VPA in raw data: varPart_raw.rdata in 5_countPeaks/variance folder                           ||
## -  Normalised counts data: ruv3_norm_prom.rdata in 5_countPeaks/normCounts folder               ||
## -  VPA in normalised counts data: varPart_ruv3_norm_prom.rdata in 5_countPeaks/variance folder  ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R > 4.3, R libraries: variancePartition, DESeq2, ruvIIInb,SummarizedExperiment,scater         ||
## - A csv file in the metadata folder with samples used for group peak calling:                   || 
##   samplesGroupAnalysis.csv. This is produced in STEP 7.0                                        ||
## - A csv file in the metadata folder with samples that passed stage 1 and 2 of QC:               || 
##   passS1S2Status.csv. This is produced in STEP 6.3                                              ||
## - (Optional) A txt file in the metadata folder with samples that need to be left out:           || 
##   leaveOutSamples.txt. This has to be manually made                                             ||
## ================================================================================================##

## ==========##
##   SET-UP  ##
## ==========##

args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
setPeaks<-args[1]
source(configFile)

library(variancePartition)
library(DESeq2)
library(ruvIIInb)
library(SummarizedExperiment)
library(scater)

## ============= ##
##   FUNCTIONS   ##
## ============= ##

## Load counts in peaks files and match to pheno samples data
loadCounts <- function(counts){
  colnames(counts) <- lapply(colnames(counts),function(x) gsub("*.filt.nodup.bam", "",x))
  counts <- counts[, colnames(counts) %in% pheno$sampleID]
  counts <- counts[,match(pheno$sampleID,colnames(counts))]
  counts <- counts[, colSums(is.na(counts)) == 0]
  return(counts)
}

## ===============##
##   GATHER DATA  ##
## ===============##

setPeaks <- if (setPeaks == "PROM") "prom" else "all"

## Load metadata and information about samples
metadata<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)

samplesPassed <- read.csv(file.path(paste0(metaDir,"/passS1S2Status.csv")), sep = ",",stringsAsFactors = FALSE)
samplesPassed <- samplesPassed[samplesPassed$QCS1 == TRUE & samplesPassed$QCS2 == TRUE,]
qcStats <- read.csv(file.path(paste0(metaDir,"/samplesGroupAnalysisQC1.csv")), sep = ",",stringsAsFactors = FALSE)


## If there are samples that we want to leave out at this stage, we check for this file if exists and adjusts the samples to be used.
if(file.exists(paste0(metaDir,"/leaveOutSamples.txt"))){
  loSamples <- read.table(file.path(paste0(metaDir,"/leaveOutSamples.txt")))[,1]
  samplesStatus <- samplesPassed[!samplesPassed$sampleID %in% loSamples,]
  pheno <- metadata[metadata$sampleID %in% samplesStatus$sampleID,]
  samplesStatus <- samplesStatus[samplesStatus$sampleID %in% pheno$sampleID,]
  qcStats <- qcStats[qcStats$Sample %in% samplesStatus$sampleID,]
}else{
  pheno <- metadata[metadata$sampleID %in% samplesPassed$sampleID,]
  qcStats <- qcStats[qcStats$Sample %in% samplesPassed$sampleID,]

}

## ============== ##
##   LOAD COUNTS  ##
## ============== ##
if(file.exists(paste0(countsDir,"/Counts/peakCounts_",setPeaks,".rdata"))) {
  print("Counts in peaks found.")
  load(paste0(countsDir,"/Counts/peakCounts_",setPeaks,".rdata"))
  counts <- loadCounts(fc_ctPeaks$count)
  pheno <- pheno[pheno$sampleID %in% colnames(counts),]
  colnames(counts) <- pheno$sampleCode
  qcStats <- qcStats[qcStats$Sample %in% pheno$sampleID,]
} else {
  stop("Counts in peaks not found. Please go back to STEP 7.3 COUNTS to obtain these.")
}


## ================================= ##
##   VARIANCE PARTITION RAW COUNTS   ##
## ================================= ##

## Add braak status as a categorical variable
braak <-c()
for(i in 1:nrow(pheno)){
  if(as.numeric(pheno$clinical[i] < 5)){
    braak<- c(braak,"low-braak")
  }else{
    braak<- c(braak,"high-braak")
  }
}
pheno$braak <- braak
colData <- data.frame(Fraction=as.factor(pheno$fraction), ID=pheno$sampleCode, Cohort=as.factor(pheno$cohort), 
                       Age=pheno$age, Batch = as.factor(pheno$sequencingBatch), 
                       Braak =as.factor(pheno$braak),Gender =as.factor(pheno$gender),
                       NRF=qcStats$NRF, NFr = qcStats$NF, MNC = qcStats$Mono, MULTIMOD = qcStats$p.value, ALIGNED_READS=qcStats$overall_alignment_rate/100)

rownames(colData) <- colData$ID

# Define formula
form <- ~ (1|Fraction) + (1|Batch)+(1|Cohort)+Age+ NRF+ NFr + MULTIMOD + ALIGNED_READS + (1|Gender)

# Run variancePartition analysis in raw data
varPart.raw <- fitExtractVarPartModel(counts, form, colData)
save(varPart.raw, file = paste0(countsDir,"/variance/varPart_raw",setPeaks,".rdata"))

if(file.exists(paste0(countsDir,"/variance/varPart_raw",setPeaks,".rdata"))) print("Variance Partition Analysis of raw data done")

## ==================== ##
##   NORMALISE COUNTS   ##
## ==================== ##

print("Counts will be normalised using RUV-III-NB method")
print(paste0("k parameter chosen is: ",as.character(k)))

## Create counts object 
se0 <- SummarizedExperiment(assays=list(counts=counts),
                            colData=colData)
se0 <-addPerFeatureQCMetrics(x = se0)

## Filter peaks with low counts
se0<-subset(se0, rowData(se0)$mean>0 )

dds <- DESeqDataSetFromMatrix(countData = assays(se0)$counts,
                              colData = colData(se0),
                              design = ~ Fraction)
dds <- DESeq(dds)
res <- results(dds)
res_sorted <- res[order(res$padj, decreasing = TRUE),]
least_significant_peaks <- rownames(head(res_sorted, 1000))

rowData(se0)$ctlLogical<-rownames(assays(se0)$counts) %in% least_significant_peaks
assayNames(se0) <- c("counts")
M <- matrix(0,ncol(assays(se0)$counts),length(unique(se0$Fraction)))
rownames(M) <- colnames(assays(se0)$counts)
colnames(M) <- unique(se0$Fraction)
for(CL in colData(se0)$Fraction){
  M[which(se0$Fraction==CL),CL] <- 1
}

ruviii.norm <-fastruvIII.nb(Y=assays(se0)$counts,M=M,ctl=rowData(se0)$ctlLogical,k=k, batch=as.numeric(colData(se0)$Batch))
save(ruviii.norm, file = paste0(countsDir,"/normCounts/ruv3_norm_",setPeaks,".rdata"))
if(file.exists(paste0(countsDir,"/normCounts/ruv3_norm_",setPeaks,".rdata"))) print(paste0("Count data has been normalised using RUV-III-NB with k=",k))

## ================================= ##
##   VARIANCE PARTITION RAW COUNTS   ##
## ================================= ##

logPAC<- assay(ruviii.norm, "logPAC")
logPAC <- logPAC[,colnames(logPAC) %in% pheno$sampleCode]
logPAC <- as.matrix(logPAC)
logPAC <- logPAC[,match(pheno$sampleCode,colnames(logPAC))]

varPart.ruv <- fitExtractVarPartModel(as.matrix(logPAC), form, colData)
save(varPart.ruv,  file = paste0(countsDir,"/variance/varPart_ruv3_norm_",setPeaks,".rdata"))
if(file.exists(paste0(countsDir,"/variance/varPart_ruv3_norm_",setPeaks,".rdata"))) print(paste0("Variance Partition analysis of normalised counts using RUV-III-NB with k=",k, " is done"))

