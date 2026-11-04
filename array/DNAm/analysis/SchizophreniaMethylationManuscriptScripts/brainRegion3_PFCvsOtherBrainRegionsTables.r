##---------------------------------------------------------------------#
##
## Title: PFC vs Other Brain Region Tables
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(dplyr)

HandScellType <- "NeuN" # "NeuN", "SOX10" are the only ones in all 3 brain regions
PFCcellType <- "NeuN+"

dataDir <- "Schizophrenia/brainRegionData/"
resDir <- file.path(dataDir, "4_analysis/EWAS")
outPath <- file.path(dataDir, "4_analysis/plots")

v2manifest <- "references/EPICArray/EPICv2_reannotated_manifest_v1.0.csv"

hippPath <- file.path(resDir, paste0(HandScellType, "HippocampusLM.rdata"))
striatumPath <- file.path(resDir, paste0(HandScellType,"StriatumLM.rdata"))
newPFCPath <- file.path(resDir, "NeuN+PFCLM.rdata")

PFCpath <- file.path(dataDir, "/4_analysis/EWAS/new/allCellType_LM_FDR.rdata")


#----------------------------------------------------------------------#
# LOAD AND PREPARE Striatum and Hippocampus Data
#----------------------------------------------------------------------#

# load epicv2 manifest to match IlmnID to EPICv1_Loci
manifest <- read.csv(v2manifest, stringsAsFactors = FALSE)
manifest <- manifest[, c("IlmnID", "EPICv1_Loci")]
colnames(manifest)[2] <- "ProbeID"

# load hippocampus and striatum results and combine
load(hippPath)
hipp <- as.data.frame(outtab)
rm(outtab)

load(striatumPath)
striatum <- as.data.frame(outtab)
rm(outtab)

regionRes <- cbind(hipp, striatum)
regionRes$IlmnID <- rownames(regionRes)

# reorder columns for combining with PFC res
matchData <- regionRes[, c(3,1,2,6,4,5,7)]
colnames(matchData) <- c("P:Hipp", "Coeff:Hipp", "SE:Hipp",
                         "P:Striatum", "Coeff:Striatum", "SE:Striatum",
                         "IlmnID")

# join with manifest to get EPICv1_Loci
matchData <- dplyr::left_join(matchData, manifest, by = "IlmnID")
matchData <- matchData[!duplicated(matchData$ProbeID),] # remove duplicates

# load new PFC EWAS results
load(newPFCPath)
newPFC <- as.data.frame(outtab)
rm(outtab)

newPFC$ProbeID <- rownames(newPFC)
newPFC <- newPFC[, c("ProbeID", "NullModel_SCZ_P", "NullModel_SCZ_coeff", "NullModel_SCZ_SE")]
colnames(newPFC) <- c("ProbeID", "P:matchedPFC", "Coeff:matchedPFC", "SE:matchedPFC")


#----------------------------------------------------------------------#
# LOAD AND PREPARE PFC Data
#----------------------------------------------------------------------#

# load PFC results files - annotated with gene info and snps/xhyb probes filtered out
load(PFCpath)

# remove 30 probes which have '*' chr and start values of 0/-1
for(i in 1:length(res)){
  res[[i]] <- res[[i]][res[[i]]$chrm != "*", ]   
}

cellRes <- res[[PFCcellType]]


#----------------------------------------------------------------------#
# Join and save
#----------------------------------------------------------------------#

cellRes$ProbeID <- rownames(cellRes)
fullRes <- dplyr::left_join(cellRes, matchData, by = "ProbeID")
fullRes <- dplyr::left_join(fullRes, newPFC, by = "ProbeID")

# subset to only useful columns
colsToKeep <- c("ProbeID","NullModel_SCZ_coeff", "NullModel_SCZ_SE", "NullModel_SCZ_P",
                "Coeff:Hipp", "SE:Hipp", "P:Hipp",                    
                "Coeff:Striatum", "SE:Striatum", "P:Striatum",
                "Coeff:matchedPFC", "SE:matchedPFC", "P:matchedPFC",
                "chrm", "start", "GeneNames")

fullRes <- fullRes[, colsToKeep]

colnames(fullRes)[2:4] <- c("Coeff:fullPFC", "SE:fullPFC", "P:fullPFC")




# save all probes
write.csv(fullRes, file = file.path(outPath, paste0(PFCcellType,"FullResultsAllBrainRegionsAllProbes.csv")), row.names = FALSE)

# save discovery DMPs
discRes <- fullRes[fullRes$`P:fullPFC` < 1e-04, ]
write.csv(discRes, file = file.path(outPath, paste0(PFCcellType,"FullResultsAllBrainRegionsAllProbesDiscoveryThreshold.csv")), row.names = FALSE)

# save with NAs removed
fullRes <- fullRes[complete.cases(fullRes), ] # remove rows with NAs
write.csv(fullRes, file = file.path(outPath, paste0(PFCcellType,"FullResultsAllBrainRegions.csv")), row.names = FALSE)


