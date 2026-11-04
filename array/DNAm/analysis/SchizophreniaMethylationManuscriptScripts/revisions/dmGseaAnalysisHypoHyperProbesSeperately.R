##---------------------------------------------------------------------#
##
## Title:dmGsea
##
##
## Input files: EWAS results file (neun only, split into hyper and hypo dmps sets)
##              betas matrix
##              
##
## Output files: results tables for gsGene, GO_BP and KEGG
##
## Next script: dmGseafilterRes.R
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



dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resData <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")

outDir <- file.path(dataDir, "4_analysis/dmGseaHypoHyperAnalysis/")


source(file.path(dataDir, "scripts/dmGsea_gsGene.R"))
source(file.path(dataDir, "scripts/dmGsea_utils.R"))
source(file.path(dataDir, "scripts/dmGsea_gsRank.R"))
source(file.path(dataDir, "scripts/dmGsea_gsProbe.R"))
source(file.path(dataDir, "scripts/dmGsea_gsPG.R"))


cellTypes <- "NeuN+"

#----------------------------------------------------------------------#
# LOAD AND FORMAT DATA
#----------------------------------------------------------------------#

load(resData)
neun <- res[[cellTypes]]
neun$Name <- row.names(neun)


# get beta values
load(normData)

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to cell types used
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]

# subset betas matrix to same samples
celltypeNormbeta <- celltypeNormbeta[, QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]



# create hypo objects
hypo.p <- neun %>% 
  dplyr::filter(NullModel_SCZ_coeff < 0) %>%
  dplyr::select(Name, NullModel_SCZ_P) %>%
  dplyr::rename(p = NullModel_SCZ_P)
 
hypoNormbeta <- celltypeNormbeta[rownames(hypo.p),]



# create hyper objects
hyper.p <- neun %>% 
  dplyr::filter(NullModel_SCZ_coeff > 0) %>%
  dplyr::select(Name, NullModel_SCZ_P) %>%
  dplyr::rename(p = NullModel_SCZ_P)

hyperNormbeta <- celltypeNormbeta[rownames(hyper.p),]



#----------------------------------------------------------------------#
# Run gsGene()
#----------------------------------------------------------------------#

setwd(file.path(outDir, "gsGene"))


#===================== Hypomethylated probes =====================#
#=================================================================#


# Threshold method KEGG
gsGeneLocal(probe.p = hypo.p, 
            Data4Cor = hypoNormbeta, 
            arrayType = "EPIC", 
            ncore = 2,
            method = "Threshold",
            gSetName = "KEGG",
            outfile = "hypo_gsGeneThreshold",
            outGenep = FALSE)



# Threshold method goBP
gsGeneLocal(probe.p = hypo.p, 
            Data4Cor = hypoNormbeta, 
            arrayType = "EPIC", 
            ncore = 2,
            method = "Threshold",
            gSetName = "GO|BP",
            outfile = "hypo_gsGeneThreshold",
            outGenep = FALSE)



#===================== Hypermethylated probes =====================#
#=================================================================#


# Threshold method KEGG
gsGeneLocal(probe.p = hyper.p, 
            Data4Cor = hyperNormbeta, 
            arrayType = "EPIC", 
            ncore = 2,
            method = "Threshold",
            gSetName = "KEGG",
            outfile = "hyper_gsGeneThreshold",
            outGenep = FALSE)



# Threshold method goBP
gsGeneLocal(probe.p = hyper.p, 
            Data4Cor = hyperNormbeta, 
            arrayType ="EPIC", 
            ncore = 2,
            method = "Threshold",
            gSetName = "GO|BP",
            outfile = "hyper_gsGeneThreshold",
            outGenep = FALSE)
  
  

 