## summary figure of peaks call per fraction

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(data.table)
library(RColorBrewer)
library(ChIPseeker)
library(clusterProfiler)
library("nucleR")


annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature="gene")

## load called peaks

neun.bed<-"NEUN_peaks.broadPeak"
sox.bed<-"SOX_peaks.broadPeak"
tot.bed<-"TOTAL_peaks.broadPeak"
neun.peaks <- toGRanges(neun.bed, format="BED", header=FALSE, skip=1)
sox.peaks <- toGRanges(sox.bed, format="BED", header=FALSE, skip=1)
tot.peaks <- toGRanges(tot.bed, format="BED", header=FALSE, skip=1)

## plot ChIP peaks coverage
pdf("Plots/PeakChrDistribution.pdf")
covplot(neun.peaks, weightCol="score")
covplot(sox.peaks, weightCol="score")
covplot(tot.peaks, weightCol="score")
dev.off()


neun.peaks.anno <- annotatePeak(neun.peaks, tssRegion=c(-3000, 3000),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
sox.peaks.anno <- annotatePeak(sox.peaks, tssRegion=c(-3000, 3000),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
tot.peaks.anno <- annotatePeak(tot.peaks, tssRegion=c(-3000, 3000),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
## write to file for use in other programs
export.bed(neun.peaks.anno@anno, filepath = "AnnotatedPeaks/NEUN_peaks_hg38", splitByChrom = FALSE, name = "NeuN+veSpecificDMRs")
export.bed(sox.peaks.anno@anno, filepath = "AnnotatedPeaks/SOX10_peaks_hg38", splitByChrom = FALSE, name = "SOX10+veSpecificDMRs")
export.bed(tot.peaks.anno@anno, filepath = "AnnotatedPeaks/TOTAL_peaks_hg38", splitByChrom = FALSE, name = "TOTALSpecificDMRs")

peakAnnoList<-list("NeuN +ve" = neun.peaks.anno, "Sox10 +ve" = sox.peaks.anno, "Total" = tot.peaks.anno)

pdf("Plots/GenomicAnnotationDistribution.pdf")
plotAnnoPie(neun.peaks.anno, main = "NeuN +ve")
plotAnnoPie(sox.peaks.anno, main = "Sox10 +ve")
plotAnnoPie(tot.peaks.anno, main = "Total")
plotAnnoBar(peakAnnoList)
dev.off()

## Visualize distribution of TF-binding loci relative to TSS
pdf("Plots/DistancetoTSS.pdf")
plotDistToTSS(peakAnnoList,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

## Functional Enrichment

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
## write significant results to 
write.csv(compKEGG@compareClusterResult, "Tables/KEGGPathwayAnalysisCalledPeaksPerFraction.csv")

dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
						 OrgDb         = org.Hs.eg.db,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
write.csv(compGO@compareClusterResult, "Tables/GOPathwayAnalysisCalledPeaksPerFraction.csv")

dotplot(compGO, showCategory = 15, title = "GO Pathway Enrichment Analysis")

## overlap of gene annotations across peaks
pdf("Plots/VennOverlapGenePeakAnno.pdf")
vennplot(genes)
dev.off()
