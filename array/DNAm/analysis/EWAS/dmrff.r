# ---------------------------------------------------------------------#
##
## Title: Regional analysis for case control within each cell type
##
## Purpose of script: identify regions of differential  DNA methylation for
## schizophrenia vs controls for each cell type separately using dmrff
##
## ---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(dmrff)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
cellType <- args[2]
refPath <- args[3]

maxgap <- 1000
p.cutoff <- 0.05

normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath <- file.path(dataDir, "4_analysis/EWAS")


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

load(normData)

print(paste0("running dmrff on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

load(file.path(resPath, paste0(cellType, "LM.rdata")))


#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#

crosshyb <- read.table(file.path(refPath, "CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter <- read.csv(file.path(refPath, "EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes <- read.table(file.path(refPath, "SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
crosshyb2 <- read.csv(file.path(refPath, "Pidsley_SM1.csv"), stringsAsFactors = FALSE)
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
remove <- intersect(remove, rownames(outtab))

outtab <- outtab[-match(remove, rownames(outtab)), ]


# add gene annotation
probeAnnot <- read.table(file.path(refPath, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
probeAnnot <- probeAnnot[match(rownames(outtab), probeAnnot$probeID), ]
probeAnnot$chrm <- gsub("chr", "", probeAnnot$chrm)
probeAnnot$start <- probeAnnot$start+1

outtab<- cbind(outtab, probeAnnot[, c("chrm", "start", "GeneNames", "GeneClasses")])
outtab<- outtab[which(outtab[,"chrm"] != "Y"),]
outtab<- outtab[which(outtab[,"chrm"] != "*"),]


celltypeNormbeta<-celltypeNormbeta[rownames(outtab),]


#----------------------------------------------------------------------#
# RUN DMR ANALYSIS
#----------------------------------------------------------------------#

## find candidate regions
candidates <- dmrff.candidates(
        estimate=outtab$FullModel_SCZ_coeff,
        p.value=outtab$FullModel_SCZ_P,
        chr=outtab$chrm, 
        pos=outtab$start,
        maxgap=maxgap,
        p.cutoff=p.cutoff,
        verbose=T)

write.csv(candidates, file = file.path(resPath, "Tables", paste0(cellType,"CollapsedRegions_0.05.csv")))

candidates <- dmrff.candidates(
        estimate=outtab$FullModel_SCZ_coeff,
        p.value=outtab$FullModel_SCZ_P,
        chr=outtab$chrm, 
        pos=outtab$start,
        maxgap=maxgap,
        p.cutoff=1e-6,
        verbose=T)

write.csv(candidates, file = file.path(resPath, "Tables", paste0(cellType,"CollapsedRegions_1e-6.csv")))


dmrs <- dmrff(estimate = outtab$FullModel_SCZ_coeff, 
    se=outtab$FullModel_SCZ_SE, 
    p.value = outtab$FullModel_SCZ_P,
    methylation = celltypeNormbeta,
    chr = outtab$chrm,
    pos = outtab$start,
    maxgap = maxgap,
    verbose = T, 
    p.cutoff = p.cutoff)

write.csv(dmrs, file = file.path(resPath, "Tables", paste0(cellType,"dmrff_0.05.csv")))