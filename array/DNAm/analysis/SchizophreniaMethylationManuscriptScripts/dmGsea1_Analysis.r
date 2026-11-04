##---------------------------------------------------------------------#
##
## Title:dmGsea
##
##
## Input files: EWAS results file (neun only)
##              betas matrix
##              
##
## Output files: results tables for gsGene and gsProbe
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



dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resData <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")


outDir <- file.path(dataDir, "4_analysis/dmGseaV3/")

# load edited package scripts (bug in original)
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
            gSetName="KEGG",
            outfile="gsGeneThreshold",
            outGenep=TRUE)



# Threshold method goBP
gsGeneLocal(probe.p = neun.p, 
            Data4Cor=celltypeNormbeta, 
            arrayType="EPIC", 
            ncore=2,
            method="Threshold",
            gSetName="GO|BP",
            outfile="gsGeneThreshold",
            outGenep=TRUE)
