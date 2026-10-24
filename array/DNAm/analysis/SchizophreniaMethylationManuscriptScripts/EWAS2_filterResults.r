##---------------------------------------------------------------------#
##
## Title: Filter results to exclude SNPs and X-hyb probes
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- "Schizophrenia"
refPath <- "references/EPICArray"

normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resPath <- file.path(dataDir, "4_analysis/EWAS/new")

cellTypes <- c("Double-", "NeuN+", "Sox10+", "Total")


#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#

setwd(dataDir)


# combine EWAS output into a list object
EWASresults<-c(list.files("4_analysis/EWAS/new", pattern = "LM.rdata", full.names = TRUE))

res<-list()
for (each in EWASresults) {
    load(each)
    res[[each]] <- outtab
}
rm(outtab)

names(res) <- cellTypes



#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#

crosshyb <- read.table(file.path(refPath, "CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter <- read.csv(file.path(refPath, "EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes <- read.table(file.path(refPath, "SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
snpProbes2 <- read.csv(file.path(refPath, "Pidsley_SM4.csv"), stringsAsFactors = FALSE)
snpProbes3 <- read.csv(file.path(refPath, "Pidsley_SM5.csv"), stringsAsFactors = FALSE)
snpProbes4 <- read.csv(file.path(refPath, "Pidsley_SM6.csv"), stringsAsFactors = FALSE)

snpProbes <- snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01), ]
snpProbes2 <- snpProbes2[which(snpProbes2$AF > 0.01), ]
snpProbes3 <- snpProbes3[which(snpProbes3$AF > 0.01), ]
snpProbes4 <- snpProbes4[which(snpProbes4$AF > 0.01), ]

dist <- cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist <- apply(dist, 1, min)
snpProbes4 <- snpProbes4[which(dist <= 10), ]

remove <- unique(c(tofilter$IlmnID, crosshyb[, 1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))
remove <- intersect(remove, rownames(res[[1]]))

if (length(remove) > 0) {
    for (i in 1:4) {
        res[[i]] <- res[[i]][-match(remove, rownames(res[[i]])), ]
    }
}


# add gene annotation
probeAnnot <- read.table(file.path(refPath, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
probeAnnot <- probeAnnot[match(rownames(res[[1]]), probeAnnot$probeID), ]
colnames(probeAnnot)[which(colnames(probeAnnot) == "CHR")] <- "chrm"
probeAnnot[["chrm"]] <- gsub("chr", "", probeAnnot[["chrm"]], fixed = TRUE)
probeAnnot$start <- probeAnnot$start+1

for(i in 1:4){
    res[[i]]<- cbind(res[[i]], probeAnnot[, c("chrm", "start", "GeneNames", "GeneClasses", "CGI", "CGIPosition")])
    res[[i]]<- res[[i]][which(res[[i]][,"chrm"] != "Y"),]
}


#----------------------------------------------------------------------#
# ADD IN FDR CORRECTION COLUMNS
#----------------------------------------------------------------------#

for(i in 1:length(EWASresults)){

    # add in FDR column
    res[[i]] <- as.data.frame(res[[i]])
    res[[i]][,"FDR_NullModel"] <- p.adjust(res[[i]][,"NullModel_SCZ_P"], method = "fdr")
    res[[i]][,"FDR_FullModel"] <- p.adjust(res[[i]][,"FullModel_SCZ_P"], method = "fdr")

}

#----------------------------------------------------------------------#
# SAVE
#----------------------------------------------------------------------#

save(res, file = file.path(dataDir, "4_analysis/EWAS/new/allCellType_LM_FDR.rdata"))
neunOnly<-res[[2]]
save(neunOnly, file = file.path(dataDir, "4_analysis/EWAS/new/NeuN+LM_FDR.rdata"))

