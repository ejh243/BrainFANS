#---------------------------------------------------------------------#
##
## Title: Regional analysis for case control within each cell type
##
## Purpose of script: identify regions of differential  DNA methylation for
## schizophrenia vs controls for each cell type separately using DMRcate
##
## Note this script uses a custom version of DMRplot that can't be added to 
## this repository due to license restrictions
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(DMRcate)

#----------------------------------------------------------------------#
# DEFINE FUNCTIONS
#----------------------------------------------------------------------#

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
cellType <- args[2]
refPath <- args[3]

normData<-file.path(dataDir, "3_normalised/normalised.rdata")
cellCompData<-file.path(dataDir, "4_analysis", "EstimatedNeuralCellComposition.rdata")
resPath<-file.path(dataDir, "4_analysis/EWAS")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)
load(cellCompData)

print(paste0("running DMRcate on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]


# add cell composition proportions & CETYGO scores

if (cellType == "Double-"){
   QCmetrics$Cell.Proportion <- predPropBest[QCmetrics$Basename, "NeuNNeg_Sox10Neg_IRF8Pos"]
} else if (cellType == "NeuN+"){
  QCmetrics$Cell.Proportion <- predPropBest[QCmetrics$Basename, "NeuNPos_SOX6Pos"]
} 

QCmetrics$CETYGO <- predPropBest[QCmetrics$Basename, "CETYGO"]

cols<-c("Basename", "Phenotype", "Cell.Proportion", "CETYGO", "CCDNAmAge", "Sex", "Tissue.Centre")
cols <- cols[cols %in% colnames(QCmetrics)]
QCmetrics<-na.omit(QCmetrics[, cols])

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#

crosshyb<-read.table(file.path(refPath, "CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter<-read.csv(file.path(refPath, "EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes<-read.table(file.path(refPath, "SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
crosshyb2<-read.csv(file.path(refPath, "Pidsley_SM1.csv"), stringsAsFactors = FALSE)
snpProbes2<-read.csv(file.path(refPath, "Pidsley_SM4.csv"), stringsAsFactors = FALSE)
snpProbes3<-read.csv(file.path(refPath, "Pidsley_SM5.csv"), stringsAsFactors = FALSE)
snpProbes4<-read.csv(file.path(refPath, "Pidsley_SM6.csv"), stringsAsFactors = FALSE)

snpProbes<-snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01),]
snpProbes2<-snpProbes2[which(snpProbes2$AF > 0.01),]
snpProbes3<-snpProbes3[which(snpProbes3$AF > 0.01),]
snpProbes4<-snpProbes4[which(snpProbes4$AF > 0.01),]

dist<-cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist<-apply(dist, 1, min)
snpProbes4<-snpProbes4[which(dist <=10),]

remove<-unique(c(tofilter$IlmnID, crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))
remove<-intersect(remove, rownames(celltypeNormbeta))

celltypeNormbeta <- celltypeNormbeta[-match(remove, rownames(celltypeNormbeta)),]


#----------------------------------------------------------------------#
# DEFINE ANALYSIS
#----------------------------------------------------------------------#

if("Cell.Proportion" %in% cols){
  design <- model.matrix(~ Phenotype + Cell.Proportion + CETYGO + CCDNAmAge + Sex + Tissue.Centre, QCmetrics)
} else {
 design <- model.matrix(~ Phenotype + CETYGO + CCDNAmAge + Sex + Tissue.Centre, QCmetrics)
}
siteAnnotation <- cpg.annotate("array", celltypeNormbeta, what = "B", arraytype = "EPIC", 
analysis.type = "differential", design = design, coef = 2)

# compare to our LM results
load(file.path(resPath, paste0(cellType, "LM.rdata")))
outtab<-outtab[names(siteAnnotation@ranges),]

pdf(file.path(resPath, "Plots", paste0("ScatterplotLMDMRCate", cellType, ".pdf")), width = 10, height = 5)
par(mfrow = c(1,2))
plot(-log10(siteAnnotation@ranges$ind.fdr), -log10(outtab[,"FullModel_SCZ_P"]), 
  xlab = "DMRcate -log10P", ylab = "LM -log10P", pch = 16, 
  col = gg_color_hue(2)[as.factor(siteAnnotation@ranges$is.sig)])
plot(siteAnnotation@ranges$diff, outtab[,"FullModel_SCZ_coeff"], 
  xlab = "DMRcate mean diff", ylab = "LM Mean diff", pch = 16, 
  col = gg_color_hue(2)[as.factor(siteAnnotation@ranges$is.sig)])
abline(v = 0)
abline(h = 0)
dev.off()

dmrcoutput <- dmrcate(siteAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")

cols<-as.factor(QCmetrics$Phenotype)
names(cols)<-QCmetrics$Basename

pdf(file.path(resPath, "Plots", paste0("DMRcate_DMRs_", cellType, ".pdf")))
for(i in 1:length(results.ranges)){
DMR.plot (ranges = results.ranges, 
                     dmr = i, 
                     CpGs = celltypeNormbeta, 
                     phen.col = cols,
                     genome = "hg19",
                     labels = names(ranges),
                     flank = 50,
                     boxplot.flank = min(c(5, round(width(results.ranges)[i]*0.05))),
                     extra.ranges = NULL) 
}
dev.off()

save(dmrcoutput, file = file.path(resPath, paste0(cellType,"DMRCate.rdata")))
write.csv(results.ranges, file.path(resPath, "Tables", paste0("DMRcate_DMRs_", cellType, ".csv")))



