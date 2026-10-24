## ---------------------------------------------------------------------#
##
## Title: CombP
##
## Purpose of script: Run compP to identify DMRs
##
## ---------------------------------------------------------------------#


# The input should be a data frame with column names "chr","start", "end","p", and "probe", 
# indicating chromosome number, start position, end position, P value and probe name

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ENmix) # contains the combp function


#----------------------------------------------------------------------#
# Set up
#----------------------------------------------------------------------#

dataDir <- "Schizophrenia"
outPath <- file.path(dataDir, "4_analysis/DMRs")
normData <- file.path(dataDir, "/3_normalised/normalised.rdata")

load(file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata"))
load(normData)

setwd(outPath)


#----------------------------------------------------------------------#
# Prepare input data
#----------------------------------------------------------------------#

# subset betas matrix and pheno to only include samples in neun results

#subset to neun samples only
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == "NeuN+"),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename] # 161 samples


# get neun EWAS results
neun <- res[[2]]

# set up cols as needed for combp
neun$probe <- rownames(neun)
neun <- neun[, c("chrm", "start", "start", "NullModel_SCZ_P", "probe")]
colnames(neun) <- c("chr", "start", "end", "p", "probe")


#----------------------------------------------------------------------#
# Run combp
#----------------------------------------------------------------------#

combp(neun,dist.cutoff=1000,bin.size=310,seed=0.18,
      region_plot=F,mht_plot=F,nCores=10,verbose=TRUE)



#----------------------------------------------------------------------#
# Define function to annotate DMRs
#----------------------------------------------------------------------#


library(pbapply) # for progress bar

# inputs:

# EWASres is EWAS results table annotated with gene names, chr and mapinfo
# dmrs is combp output data frame
# EScol is the name of the column in EWASres containing effect sizes
# Pcol is the name of the column in EWASres containing p-values


# function to collapse multiple gene names into one value
uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }


# function to annotate combp DMRs prior to plotting 
annot_dmrs <- function(DMRrow, EWASres, EScol, Pcol){
  
  ewasCols <- c("CHR", "MAPINFO", "Gene")
  if(!all(ewasCols %in% colnames(EWASres))){
    stop(paste("EWAS results must contain the following columns:", paste(ewasCols, collapse = ", ")))
  }
  if(!(EScol %in% colnames(EWASres))){    
    stop(paste("Effect size column", EScol, "not found in EWAS results"))
  }
  if(!(Pcol %in% colnames(EWASres))){   
    stop(paste("P-value column", Pcol, "not found in EWAS results"))
  }
  
  
  x.chr <- as.character(DMRrow[colnames(dmrs)=='chr'])
  x.start <- as.numeric(DMRrow[colnames(dmrs)=='start'])
  x.end <- as.numeric(DMRrow[colnames(dmrs)=='end'])
  dmr_stats <- EWASres[as.character(EWASres$CHR)==x.chr & EWASres$MAPINFO>=x.start & EWASres$MAPINFO<=x.end,]
  max.DMP.Estimate <- dmr_stats[which.max(abs(dmr_stats[, EScol])), EScol]
  min.DMP.P <- min(dmr_stats[,Pcol])
  Gene <- paste(unique(dmr_stats$Gene[dmr_stats$Gene!='']),collapse=';') # collapse multiple genes into one value
  Gene <- unlist(lapply(Gene, uniqueAnno))
  return(c(max.DMP.Estimate, min.DMP.P, Gene))
}



# function that wraps it all together

annotateDMRs <- function(dmrs, EWASres, EScol, Pcol){
  
  print(paste0("Annotating ", nrow(dmrs), " DMRs"))
  
  annots <- t(pbapply(dmrs, 1, annot_dmrs, EWASres=EWASres, EScol=EScol, Pcol=Pcol)) # add progress bar
  
  # Convert to data frame and set column names
  annots <- as.data.frame(annots, stringsAsFactors = FALSE)
  colnames(annots) <- c('max.DMP.Estimate','min.DMP.P','Gene')
  
  # Combine with original DMRs
  dmrs <- cbind(dmrs, annots)
  
  return(dmrs)
}



#----------------------------------------------------------------------#
# Annotate results
#----------------------------------------------------------------------#

# read in combp results
dmrs <- read.csv("resu_combp.csv")   

# remove regions with < 3 probes
dmrs <- dmrs[dmrs$nprobe > 2,]

# get full EWAS results table and ensure cols match with function
neun <- res[[2]]
colnames(neun)[which(colnames(neun) == "chrm")] <- "CHR"
neun$MAPINFO <- neun$start
colnames(neun)[which(colnames(neun) == "GeneNames")] <- "Gene"

# annotate DMRs
DMRsAnnotated <- annotateDMRs(dmrs, EWASres=neun, EScol="NullModel_SCZ_coeff", Pcol="NullModel_SCZ_P")

# save results
write.csv(DMRsAnnotated, file=file.path(outPath, "combp1e04Annotated.csv"), row.names=F)

