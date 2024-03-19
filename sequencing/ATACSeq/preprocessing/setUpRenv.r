## ============================================================================================##
##                     ATAC-seq pipeline STEP 0.2: Set up R environment                        ##
## ============================================================================================##
## EXECUTION: Rscript ./sequencing/ATACSeq/preprocessing/setUpRenv.sh  <user>                  ||
## - execute from scripts directory                                                            ||
##                                                                                             ||
## DESCRIPTION: This script checks if required R packages are installed and install them if not//
##                                                                                             ||
## INPUTS:                                                                                     ||
## - <project> : project on which analysis is being run                                        ||
##                                                                                             ||
## REQUIRES:                                                                                   ||
## - R/4.2.1-foss-2022a                                                                        ||
## ============================================================================================##

args <- commandArgs()
user <-args[6]

##Vectors with all required R packages in the pipeline
packagesBase <- c("diptest", "ptest", "plyr","dplyr","ggplot2","ggpubr","RColorBrewer",  "pheatmap", "gridExtra", "FME", "data.table", "readr","scales","tibble","kableExtra","tidyverse","rstatix","conflicted","knitr","vioplot",
"corrplot","reshape2","cowplot")
packagesBiocM <- c("BiocManager","ATACseqQC","Rsubread","GenomicRanges","ChIPpeakAnno","ChIPseeker","TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db","edgeR","DiffBind","csaw","GenomicAlignments", "GenomicTools.fileHandler", 
"rtracklayer", "clusterProfiler","ReactomePA")

## As most packages are installed through BiocManager, we first installed this if not found
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(packagesBiocM, lib = paste0("/lustre/home/",user,"/R/x86_64-pc-linux-gnu-library/4.2") )

## Then we install the rest of packages if not installed already
install.packages(setdiff(packagesBase, rownames(installed.packages())),lib = paste0("/lustre/home/",user,"/R/x86_64-pc-linux-gnu-library/4.2"))  