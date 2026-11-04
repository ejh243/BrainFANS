##---------------------------------------------------------------------#
##
## Title: Mini Manhatten plots
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
library(data.table)
library(dplyr)
library(GenomicRanges)
library(Gviz)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(rtracklayer)
library(ComplexHeatmap)
library(biomaRt)



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#


dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resData <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")


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


# ewas data
ewas <- neun
ewas$chr <- ewas$chrm
ewas$pos <- ewas$start
ewas$P <- ewas$NullModel_SCZ_P
ewas$probeID <- row.names(ewas)
ewas$deltaBeta <- ewas$NullModel_SCZ_coeff



# enchancers etc.
download.file(
  "https://www.encodeproject.org/files/ENCFF420VPZ/@@download/ENCFF420VPZ.bed.gz",
  "neuron_enhancers.bed.gz"
)


enh_df <- read.table("neuron_enhancers.bed.gz", header = FALSE, sep = "\t", 
                     comment.char = "#", fill = TRUE)

colnames(enh_df) <- c("chr", "start", "end", "name", "score", "strand", 
                      "thickStart", "thickEnd", "color", "type")

enh_df <- enh_df[enh_df$type %in% c("pELS", "dELS", "PLS"), ]



#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#

pdf("MiniMan.pdf")

for(i in 1:16){

  top_probe <- ewas %>%
    arrange(NullModel_SCZ_P) %>%
    slice(i)

  probe <- top_probe$probeID

  vLarge <- c(7, 15)
  larger <- 10
  smaller <- c(9, 12, 16)

  window <- 30000

  if(i %in% larger){
    window <- window + 30000
  } else if (i %in% smaller){
    window <- window - 10000
  } else if (i %in% vLarge){
    window <- window + 90000
  }



  chr <- top_probe$chr
  centre <- top_probe$pos
  from <- centre - window
  to   <- centre + window

region <- ewas %>%
  dplyr::filter(chrm == !!chr, pos >= from, pos <= to) %>%
  mutate(
    logP = -log10(P),
    is_top = probeID == top_probe$probeID
  )

# label only strongest probes
label_n <- 8
region <- region %>%
  arrange(P) %>%
  mutate(label = ifelse(row_number() <= label_n | is_top, probeID, ""))



# EWAS data track
ewas_track <- DataTrack(
  range = GRanges(
    seqnames = chr,
    ranges = IRanges(start = region$pos, end = region$pos),
    probeID = region$probeID,
    deltaBeta = region$deltaBeta,
    label = ifelse(region$is_top, region$probeID, "")   # <-- add label column
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


# probe track
probe_track <- AnnotationTrack(
  range = GRanges(
    seqnames = chr,
    ranges   = IRanges(start = region$pos, end = region$pos),
  ),
  genome   = "GRCh38",
  chromosome = chr,
  name     = "EPIC probes",
  shape    = "box",
  stacking = "dense",
  showId   = TRUE,
  cex      = 0.55,
)


# gene track
grt_data <- getGeneRegionTrackForGviz(
  EnsDb.Hsapiens.v86,
  chromosome = chr,
  start = from,
  end = to
)


# Remove rows where gene starts with "RP11-"
grt_data <- grt_data[!grepl("^RP11-|^CTA-", grt_data$symbol), ]


rights <- c(3, 7, 10, 12, 13, 15 )

if(i %in% rights){
  labelPos <- "right"
} else {
  labelPos <- "below"
}


gene_track <- GeneRegionTrack(
  grt_data,
  genome = "GRCh38",
  chromosome = chr,
  name = "Genes",
  transcriptAnnotation = "symbol",
  collapseTranscripts = "meta",
  showId = TRUE,
  stacking = "squish",
  just.group = labelPos
)


# enhancer annotations

enh <- GRanges(
  seqnames = sub("^chr", "", enh_df$chr),
  ranges = IRanges(start = enh_df$start, end = enh_df$end),
  name = enh_df$name,
  type = enh_df$type
)


enh_region <- enh[seqnames(enh) == chr &
                    start(enh) <= to &
                    end(enh) >= from]


if(length(enh_region) > 0){
  enh_track <- AnnotationTrack(
    range = enh_region,
    genome = "GRCh38",
    chromosome = chr,
    name = "Regulatory",
    fill = ifelse(enh_region$type == "PLS",  "firebrick",
                  ifelse(enh_region$type == "pELS", "darkorange", "gold")),
    col  = ifelse(enh_region$type == "PLS",  "firebrick",
                  ifelse(enh_region$type == "pELS", "darkorange", "gold")),
    stacking = "dense"
  )
} else {
  enh_track <- AnnotationTrack(genome = "GRCh38", chromosome = chr, name = "Regulatory")
}


# create legend
lgd <- packLegend(
  Legend(labels = "Promoter (PLS)",
         legend_gp = gpar(fill = "firebrick"),
         title = ""),
  Legend(labels = "Proximal enhancer (pELS)",
         legend_gp = gpar(fill = "darkorange"),
         title = ""),
  Legend(labels = "Distal enhancer (dELS)",
         legend_gp = gpar(fill = "gold"),
         title = ""),
  direction = "horizontal"
)


legend_track <- CustomTrack(
  plottingFunction = function(GdObject, prepare, ...) {
    if(!prepare){
      draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
    }
    invisible(GdObject)
  },
  name = " ",        # blank name so no title bar text
  size = 0.2         # adjust height of the legend track
)




# genome axis track 
axis_track <- GenomeAxisTrack()



# create full plot
plotTracks(
  list(
    axis_track,
    ewas_track,
    probe_track,
    gene_track,
    enh_track,
    legend_track
  ),
  from = from,
  to = to,
  chromosome = chr,
  genome = "GRCh38"
)

}

dev.off()

