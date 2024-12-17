## ============================================================================================##
##                     ATAC-seq pipeline STEP 0.2: Set up R environment                        ##
## ============================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/Rscripts/installLibraries.r  <user>                 ||
## - execute from scripts directory                                                            ||
##                                                                                             ||
## DESCRIPTION: This script checks if required R packages are installed and install them if not||
##                                                                                             ||
## INPUTS:                                                                                     ||
## - <project> : project on which analysis is being run                                        ||
##                                                                                             ||
## REQUIRES:                                                                                   ||
## - R version > 4.3                                                                           ||
## ============================================================================================##

args <- commandArgs(trailingOnly=TRUE)
rlibrary <-args[1]

##Vectors with all required R packages in the pipeline
packagesBase <- c("diptest", "plyr","dplyr","ggplot2","ggpubr","RColorBrewer",  "pheatmap", "gridExtra", "FME", "data.table", "readr","scales","tibble","kableExtra","tidyverse","rstatix","conflicted","knitr","vioplot","corrplot","reshape2","cowplot", "devtools", "DT", "plotly", "ggExtra", "Matrix")
packagesBiocM <- c("BiocManager","ATACseqQC","Rsubread","GenomicRanges","ChIPpeakAnno","ChIPseeker","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db","edgeR","DiffBind","csaw","GenomicAlignments", "GenomicTools.fileHandler", "BSgenome.Hsapiens.UCSC.hg38",
"rtracklayer", "clusterProfiler","ReactomePA","DESeq2", "ChIPQC", "variancePartition", "SummarizedExperiment", "scater", "cluster")

## As most packages are installed through BiocManager, we first installed this if not found
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18",lib = rlibrary)

BiocManager::install(packagesBiocM, lib = rlibrary) 

## Then we install the rest of packages if not installed already
install.packages(setdiff(packagesBase, rownames(installed.packages())),lib = rlibrary)  
install.packages(paste0(dir,"/config/ptest_1.0-8.tar.gz",lib = rlibrary)
devtools::install_github("limfuxing/ruvIIInb",build_vignettes=FALSE)