##---------------------------------------------------------------------#
##
## Title:dmGsea
##
##
## Input files: EWAS results file (neun only)
##              betas matrix
##              
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# https://bioconductor.org/packages/devel/bioc/vignettes/dmGsea/inst/doc/dmGsea.html



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(dmGsea)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(dplyr)
library(stringr)


dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resData <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")


outDir <- file.path(dataDir, "4_analysis/dmGseaNewThresholds")

source(file.path(dataDir, "scripts/dmGsea_gsGene.R"))
source(file.path(dataDir, "scripts/dmGsea_utils.R"))
source(file.path(dataDir, "scripts/dmGsea_gsRank.R"))
source(file.path(dataDir, "scripts/dmGsea_gsProbe.R"))
source(file.path(dataDir, "scripts/dmGsea_gsPG.R"))



#----------------------------------------------------------------------#
# LOAD AND FORMAT DATA
#----------------------------------------------------------------------#

load(resData)
neun <- res[[2]]
neun$Name <- row.names(neun)

load(normData)

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to cell types used
cellTypes <- "NeuN+"
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]

# subset betas matrix to same samples
celltypeNormbeta <- celltypeNormbeta[, QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]
celltypeNormbeta <- celltypeNormbeta[rownames(neun),]


# create probe.p object
neun.p <- neun %>% 
  dplyr::select(Name, NullModel_SCZ_P) %>%
  dplyr::rename(p = NullModel_SCZ_P)


thresh = 0.01
#thresh = 0.1


#----------------------------------------------------------------------#
# Run gsGene()
#----------------------------------------------------------------------#

setwd(file.path(outDir, "gsGene"))



# Threshold method KEGG
gsGeneLocal(probe.p = neun.p, 
            Data4Cor=celltypeNormbeta, 
            arrayType="EPIC", 
            ncore=2,
            method="Threshold",
            FDRthre = thresh,
            gSetName="KEGG",
            outfile=paste0("gsGeneThresholdFDR", thresh), 
            outGenep=TRUE)



# Threshold method goBP
gsGeneLocal(probe.p = neun.p, 
            Data4Cor=celltypeNormbeta, 
            arrayType="EPIC", 
            ncore=2,
            method="Threshold",
            FDRthre = thresh,
            gSetName="GO|BP",
            outfile=paste0("gsGeneThresholdFDR", thresh), 
            outGenep=TRUE)


#----------------------------------------------------------------------#
# CREATE SUMMARY TABLE
#----------------------------------------------------------------------#

setwd(outDir)


# read in results files
newFiles <- list.files("gsGene", pattern = ".csv", full.names = T)
oldFiles <- list.files("../dmGseaV3/gsGene/", pattern = "gsGeneThreshold", full.names = T)
fileList <- c(newFiles, oldFiles)

resFiles<- lapply(fileList, function(x) read.csv(x, stringsAsFactors = F))
names(resFiles) <- tools::file_path_sans_ext(basename(fileList))


# function to extract info from names
parse_run_name <- function(nm) {
  fdr <- str_extract(nm, "(?<=FDR)[0-9.]+")
  fdr <- ifelse(is.na(fdr), "0.05", fdr)
  onto <- ifelse(str_detect(nm, "KEGG"), "KEGG", "GO_BP")
  list(FDR_threshold = fdr, Ontology = onto)
}


# create table
run_names <- names(resFiles)

summary_list <- lapply(run_names, function(run) {
  df <- resFiles[[run]]
  parsed <- parse_run_name(run)
  sig <- df %>% dplyr::filter(nGene > 1 & FDR < 0.05)
  
  data.frame(
    Ontology            = parsed$Ontology,
    FDRthreshold       = parsed$FDR_threshold,
    nSigPathways               = nrow(sig),
    percSigPathways             = round(100 * nrow(sig) / nrow(df), 1),
    median_nSigGene = median(sig$nSigGene, na.rm = TRUE),
    mean_percSigGene     = round(mean((sig$nSigGene / sig$nGene) * 100, na.rm = TRUE), 1)

  )
})

summary_tbl <- bind_rows(summary_list) %>%
  arrange(Ontology, as.numeric(FDRthreshold))


#----------------------------------------------------------------------#
# SAVE
#----------------------------------------------------------------------#

write.csv(summary_tbl, "thresholdsSummaryTable.csv")


