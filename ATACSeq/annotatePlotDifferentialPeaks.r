## annotate chip-seq peaks with nearest genes

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(data.table)
library(RColorBrewer)
library(ChIPseeker)


annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature="gene")

## load differentially called peaks

bed1<-"neun_vs_sox_c3.0_cond1.bed"
bed2<-"neun_vs_sox_c3.0_cond2.bed"
gr1 <- toGRanges(bed1, format="BED", header=FALSE, skip=1)
gr2 <- toGRanges(bed2, format="BED", header=FALSE, skip=1)

## flip score on second file so can merge peaks
gr2$score<-1/gr2$score
allpeaks<-c(gr1,gr2)
allpeaks<-sort(allpeaks, by = ~ score) # nb this is pointless as when you runte annotation commandslater they get resorted into chr order.

rm(gr1,gr2)

## summarise location of peaks
aCR<-assignChromosomeRegion(allpeaks, nucleotideLevel=FALSE, 
                           precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                           TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)

pdf("BarplotGenomicAnnotationNeunvsSox.pdf")
par(mar = c(10, 4,4,1))						   
barplot(aCR$percentage, las=3, main = "Genomic distribution of differential peaks", ylab = "% differential peaks")
dev.off()

## annotate peaks to promoter regions of genes	
# NB not all regions will be annotated					   
allpeaks.promotor <- annotatePeakInBatch(allpeaks,AnnotationData=annoData,output="nearestBiDirectionalPromoters",bindingRegion=c(-2000, 500))
## add gene symbol
allpeaks.promotor <- addGeneIDs(allpeaks.promotor,"org.Hs.eg.db",IDs2Add = "symbol", feature_id_type="entrez_id")

## annotate peaks to the nearest genes						   
allpeaks.nearest <- annotatePeakInBatch(allpeaks,AnnotationData=annoData,output="both", maxgap=5000)
## add gene symbol
allpeaks.nearest <- addGeneIDs(allpeaks.nearest,"org.Hs.eg.db",IDs2Add = "symbol", feature_id_type="entrez_id")

allpeaks.genes<-annotatePeak(allpeaks, tssRegion=c(-3000, 3000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
allpeaks.genes<-as.GRanges(allpeaks.genes)
geneIdMap<-addGeneIDs(allpeaks.genes$geneId,"org.Hs.eg.db",IDs2Add = "symbol", feature_id_type="entrez_id")

allpeaks.genes$symbol<- geneIdMap[match(allpeaks.genes$geneId,geneIdMap$entrez_id),"symbol"]


## plot examples using Gviz

axTrack <- GenomeAxisTrack()
plotCols<-brewer.pal(4, "Set1")
window<-3000

pdf("ExamplesDifferentialPeaksNeunSox.pdf", width = 10, height = 6)
for(i in c(1:6, tail(1:length(allpeaks.genes), n = 6))){
	chr <- as.character(seqnames(allpeaks.genes)[i])
	idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)
	## if a gene nearby include in the plot otherwise just had window around the peak
	if(abs(allpeaks.genes$distanceToTSS[i]) < window){
		start <- min(c(allpeaks.genes$geneStart[i],start(allpeaks.genes)[i]))
		end <- max(c(allpeaks.genes$geneEnd[i],end(allpeaks.genes)[i]))
		} else{
			start <- start(allpeaks.genes)[i]
			end <-end(allpeaks.genes)[i]
		}
	start<-start-window
	end<-end+window

	## double check if any additional diff peaks in region
	peaksInRegion<-findOverlaps(allpeaks.genes, GRanges(chr, IRanges(start, end))) 
	peaksInRegion<-allpeaks.genes[queryHits(peaksInRegion),]

	refGenes <- UcscTrack(genome="hg38", chromosome=chr, table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=end, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", fill="#8282d2",
	name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")
	cpgIslands <- UcscTrack(genome="hg38", chromosome=chr, track="cpgIslandExt", from=start, to=end,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

	## load bdg 
	neun.bdg<-fread(cmd=paste("grep",chr,"NEUN_treat_pileup.bdg"), header = FALSE, sep = "\t")
	neun.bdg<-neun.bdg[with(neun.bdg, ((V2 < end & V2 >= start) | (V3 <=end & V2>= start))),]
	neun.granges<-GRanges(seqnames = neun.bdg$V1, IRanges(start = neun.bdg$V2, end = neun.bdg$V3), score = neun.bdg$V4)
						 
	sox.bdg<-fread(cmd=paste("grep",chr,"SOX_treat_pileup.bdg"), header = FALSE, sep = "\t")
	sox.bdg<-sox.bdg[with(sox.bdg, ((V2 < end & V2 >= start) | (V3 <=end & V2>= start))),]
	sox.granges<-GRanges(seqnames = sox.bdg$V1, IRanges(start = sox.bdg$V2, end = sox.bdg$V3), score = sox.bdg$V4)

	y_lim<-c(0, max(c(sox.bdg$V4, neun.bdg$V4)))
	dTrack.neun <- DataTrack(range = neun.granges, genome = "hg38", type = "polygon", chromosome = chr, name = "NeuN+ve", col = plotCols[3], fill.mountain = c(plotCols[3], plotCols[3]), ylim = y_lim, baseline = 0)
	dTrack.sox <- DataTrack(range = sox.granges, genome = "hg38", type = "polygon", chromosome = chr, name = "Sox10+ve", col = plotCols[2], fill.mountain = c(plotCols[2], plotCols[2]), ylim = y_lim, baseline = 0)
	ht <- HighlightTrack(trackList=list(axTrack, refGenes,cpgIslands, dTrack.neun,dTrack.sox), start=start(peaksInRegion), end = end(peaksInRegion),chromosome=chr)					 
	plotTracks(list(idxTrack, ht), from=start, to = end,chromosome=chr, showTitle=TRUE)

}
					 
dev.off()
