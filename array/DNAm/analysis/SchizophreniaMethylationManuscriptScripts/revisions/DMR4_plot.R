## ---------------------------------------------------------------------#
##
## Title: Mini manhattan plots for ZIC1 DMRs
##
##
## ---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#


library(data.table)
library(dplyr)
library(GenomicRanges)
library(Gviz)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(rtracklayer)
library(ComplexHeatmap)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
dmrResData <- file.path(dataDir, "4_analysis/DMR/DMRcateAndCombpDMRsWithDMP.csv")
ewasResData <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")


outDir <- file.path(dataDir, "4_analysis/DMR")

setwd(dataDir)


#----------------------------------------------------------------------#
# LOAD AND FORMAT DATA
#----------------------------------------------------------------------#

dmrRes <- read.csv(dmrResData)

load(ewasResData)
neun <- res[[2]]
neun$Name <- row.names(neun)


load(normData)

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to cell types used
cellTypes <- "NeuN+"
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]


# EWAS input
ewas <- neun
ewas$pos <- ewas$start
ewas$P <- ewas$NullModel_SCZ_P
ewas$probeID <- row.names(ewas)
ewas$deltaBeta <- ewas$NullModel_SCZ_coeff



#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#


# First two DMRs (both annotated to ZIC1 gene)
dmr <- dmrRes[1:2, ]

## get gene info
zic1 <- ensembldb::genes(EnsDb.Hsapiens.v86,
                         filter = GeneNameFilter("ZIC1"))
zic1 <- zic1[seqnames(zic1) %in% c(1:22, "X", "Y")]

chr  <- as.character(seqnames(zic1))[1]
from <- min(start(zic1))
to   <- max(end(zic1))


region <- ewas %>%
  dplyr::filter(chrm == !!chr, pos >= from, pos <= to) %>%
  mutate(
    logP = -log10(P),
    is_top = probeID == probeID[which.min(NullModel_SCZ_P)]
  )

# label only strongest probes
label_n <- 8
region <- region %>%
  arrange(P) %>%
  mutate(label = ifelse(row_number() <= label_n | is_top, probeID, ""))


# EWAS track
ewas_track <- DataTrack(
  range = GRanges(
    seqnames = chr,
    ranges = IRanges(start = region$pos, end = region$pos),
    probeID = region$probeID,
    deltaBeta = region$deltaBeta
  ),
  data = region$logP,
  genome = "GRCh38",
  chromosome = chr,
  name = "-log10(P)",
  type = "p",
  cex = 0.8,
  col = ifelse(region$is_top, "red3", "grey30"),
  fill = ifelse(region$is_top, "red3", "grey30"),
  ylim = c(0, max(region$logP, na.rm = TRUE) * 1.15)
)


# line for sig threshold
sig_threshold <- -log10(1e-04)

threshold_track <- DataTrack(
  range = GRanges(
    seqnames = chr,
    ranges = IRanges(start = c(from, to), end = c(from, to))
  ),
  data = c(sig_threshold, sig_threshold),
  genome = "GRCh38",
  chromosome = chr,
  type = "l",
  col = "blue",
  lty = "dashed",
  lwd = 1,
  ylim = c(0, max(region$logP, na.rm = TRUE) * 1.15)
)


# Overlay threshold line on top of the EWAS track
ewas_overlay <- OverlayTrack(
  trackList = list(ewas_track, threshold_track),
  name = "-log10(P)"
)

# Add shading forDMRs
highlight_track <- HighlightTrack(
  trackList = list(ewas_overlay),
  start = dmr$start,
  end = dmr$end,
  chromosome = chr,
  col = "grey90",
  fill = "grey90"
)


deltabeta_track <- DataTrack(
  range = GRanges(
    seqnames = chr,
    ranges = IRanges(start = region$pos, end = region$pos)
  ),
  data = region$deltaBeta,
  genome = "GRCh38",
  chromosome = chr,
  name = "Effect Size",
  type = "h",
  baseline = 0,
  col.baseline = "grey40"
)


# gene track
grt_data <- getGeneRegionTrackForGviz(
  EnsDb.Hsapiens.v86,
  chromosome = chr,
  start = from,
  end = to
)

# Remove rows where gene starts with "RP11-" or "CTA-"
grt_data <- grt_data[!grepl("^RP11-|^CTA-", grt_data$symbol), ]


gene_track <- GeneRegionTrack(
  grt_data,
  genome = "GRCh38",
  chromosome = chr,
  name = "Genes",
  transcriptAnnotation = "symbol",
  collapseTranscripts = "meta",
  stacking = "squish",
  just.group = "above",
  fontsize.group = 12,
  fontcolor.group = "black",
  min.width = 1,
  min.height = 3
)


# genome axis track
axis_track <- GenomeAxisTrack()
ideo_track <- IdeogramTrack(genome = "hg38", chromosome = paste0("chr", chr))



# create and save

pdf(file.path(outDir, "combpAndDMRcateZIC1DMRMiniMan.pdf"), height = 12, width = 10)

plotTracks(
  list(
    ideo_track,
    axis_track,
    highlight_track,
    deltabeta_track,
    gene_track
  ),
  from = from,
  to = to,
  chromosome = chr,
  genome = "GRCh38",
  sizes = c(1, 1, 4, 1, 1),
  background.title = "grey40",   # darker title bar
  fontcolor.title = "white" 
)

dev.off()

