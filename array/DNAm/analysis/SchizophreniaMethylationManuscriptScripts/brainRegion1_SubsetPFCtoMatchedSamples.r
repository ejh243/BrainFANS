##---------------------------------------------------------------------#
##
## Title: Subset PFC normalised betas to only sample run in other brain regions
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# SUBSET
#----------------------------------------------------------------------#

# load brain region QCmetrics
load("Schizophrenia/brainRegionData/3_normalised/normalised.rdata")
BRQCmetrics <- QCmetrics
IDs <- unique(BRQCmetrics$Individual_ID)
IDs <- gsub("_", "-", IDs) # change "_" to "-" in IDs

# load PFC QCmetrics
load("Schizophrenia/3_normalised/normalised.rdata")
PFCQCmetrics <- QCmetrics

# subset PFC QCmetrics to match IDs in brain region QCmetrics
QCmetrics <- PFCQCmetrics[PFCQCmetrics$Indidivual.ID %in% IDs, ]
PFCids <- unique(QCmetrics$Indidivual.ID)

setdiff(IDs, PFCids) # 0, all IDs match
setdiff(PFCids, IDs) # 0, all IDs match

# subset PFC betas to just these samples
celltypeNormbeta <- celltypeNormbeta[, colnames(celltypeNormbeta) %in% QCmetrics$Basename]

save(QCmetrics, celltypeNormbeta, file = "Schizophrenia/brainRegionData/3_normalised/PFCnormalised.rdata")

