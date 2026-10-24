##---------------------------------------------------------------------#
##
## Title:create DMR table
##
##---------------------------------------------------------------------#



#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(dplyr)
library(stringr)


dataDir <- "Schizophrenia"
resFile <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")

cellType <- "NeuN+"


outDir <- file.path(dataDir, "4_analysis/DMR/")
outFile <- file.path(outDir, paste0(cellType, "DMRcateDMRs.rdat"))


setwd(outDir)


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

# load results files - annotated with gene info and snps/xhyb probes filtered out
load(resFile)
res <- res[[cellType]]
res <- res[res$chrm != "*",] # remove 30 probes which have '*' chr and start values of 0/-1

sig <- res[res$NullModel_SCZ_P < 1e-04,]

# dmp: character vector of CpG ids, e.g. rownames(dmp_table) or dmp$probe
dmp_ids <- rownames(sig)


#----------------------------------------------------------------------#
# ADD DMPS
#----------------------------------------------------------------------#

# load combp results
combp <- read.csv("combp/combp1e04Annotated.csv")

combp$isDMP <- vapply(
  strsplit(as.character(combp$probe), ";", fixed = TRUE),
  function(p) {
    hits <- intersect(trimws(p), dmp_ids)
    if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
  },
  character(1)
)

write.csv(combp, "combp/combpWithDMP.csv")


#----------------------------------------------------------------------#
# MANUALLY ADD DMRCATE DMPS
#----------------------------------------------------------------------#

# load DMRcate results
DMRcate <- read.csv("DMRcate/NeuN+DMRcateDMRsHG38.csv")


# DMRs added to results file

# chrm 10 (48604641 - 48605389) is in compb (48604640 - 48605389 
# chrm 19 (3982568 - 3983221) is in compb (3982567 - 3983221)
# chrm 20 (43513365 - 43514307) is in combp (43513364 - 43514405)
# chrm 3 (147405336 - 147406313) is in combp (147405335 - 147407040)
# chrm 17 (19748318 - 19748860) is in combp (19748317 - 19748944)
# chrm 7 (5495671 - 5496061) is in combp (5495089 - 5496303)
# chrm 11 (78676145 - 78676245) is in combp (78676144 - 78676524)
# chrm 3 (147423801 - 147423877) is in combp (147422330 - 147424628)
# chrm 10 (132306178 - 132306241) is in combp (132305556 - 132306525)



# subset to DMRs in both combp and DMRcate

tab <- read.csv("combp/combpWithDMPandDMRcateDMRs.csv")
tab <- tab[!is.na(tab$DMRcate_smoothedFDR),] # remove rows that are not both 

write.csv(tab, "DMRcateAndCombpDMRsWithDMP.csv")








