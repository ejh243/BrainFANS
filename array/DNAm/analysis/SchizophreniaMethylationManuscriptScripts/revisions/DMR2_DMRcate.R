##---------------------------------------------------------------------#
##
## Title: DMRcate
##
##
## Input files: normalised betas and QCmetrics
##
## Output files: csv of DMRs
##---------------------------------------------------------------------#



#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(DMRcate)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(rtracklayer)


dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resFile <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")

cellType <- "NeuN+"

outDir <- file.path(dataDir, "4_analysis/DMR/DMRcate")
outFile <- file.path(outDir, paste0(cellType, "DMRcateDMRs.rdat"))


setwd(outDir)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

# load results files - annotated with gene info and snps/xhyb probes filtered out
load(resFile)
res <- res[[cellType]]
res <- res[res$chrm != "*",] # remove 30 probes which have '*' chr and start values of 0/-1



load(normData)

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# subset beta matrix to analysis samples
celltypeNormbeta <- celltypeNormbeta[rownames(res), QCmetrics$Basename]


#----------------------------------------------------------------------#
# RUN DMRcate
#----------------------------------------------------------------------#

design <- model.matrix(~QCmetrics$Phenotype + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Tissue.Centre)

myAnnotation <- cpg.annotate(datatype = "array", arraytype = "EPICv1", celltypeNormbeta, analysis.type="differential", design=design, contrasts = FALSE, coef="QCmetrics$PhenotypeSchizophrenia", annotation = "ilm10b2.hg19")

DMRs = dmrcate(myAnnotation, lambda=1000, C=2)

results.ranges <- extractRanges(DMRs, genome = "hg19") 

dfRes <- as.data.frame(results.ranges)
write.csv(dfRes, file = paste0(outFile, ".csv"))


#----------------------------------------------------------------------#
#  LIFTOVER TO HG38
#----------------------------------------------------------------------#

# download https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain

chain <- import.chain("hg19ToHg38.over.chain")

lifted <- liftOver(results.ranges, chain) 

lifted_n <- lengths(lifted)
table(lifted_n)

lifted_gr <- unlist(lifted[lifted_n == 1])

dmrRes <- data.frame(
  chr   = sub("^chr", "", as.character(seqnames(lifted_gr))),
  start = start(lifted_gr),
  end   = end(lifted_gr),
  as.data.frame(mcols(lifted_gr))   # carries over Stouffer p, maxbetafc, etc. from extractRanges()
)

write.csv(dmrRes, "NeuN+DMRcateDMRsHG38.csv")

