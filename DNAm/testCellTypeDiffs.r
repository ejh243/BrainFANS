anovaModel<-function(row, celltype, age,sex){
	model<-lm(row ~ celltype + age + sex)
	null<-lm(row ~ age + sex)
	tmp<-summary(model)$coefficients
	regCoeffs<-c(t(tmp[c("celltypeNeuN +ve", "celltypeSox10 +ve"), c(1,4)]), anova(model,null)[2,6])
	return(regCoeffs)
}

library(ChIPseeker)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(doParallel)
library(RColorBrewer)


### test cell comp differences
source("rmdConfig.run1")
setwd(dataDir)
load("NormalisedQCd.rdata")

## double check PCs
#pca <- prcomp(t(celltypenormbeta))

## remove rs probes
celltypenormbeta<-celltypenormbeta[grep("rs", rownames(celltypenormbeta), invert = TRUE),]

## get probe annotation
annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))
newfData <- newfData[rownames(celltypenormbeta), ] # SNP probes will be missing, and be NA’d

## remove Y chr
celltypenormbeta<-celltypenormbeta[which(newfData$chr != "chrY"),]
newfData <- newfData[rownames(celltypenormbeta), ]

## remove cross-hybirdising & snp probes
crosshyb<-read.table(paste(refFiles, "/EPICArray/CrossHydridisingProbes_McCartney.txt", sep = ""), stringsAsFactors = FALSE)
snpProbes<-read.table(paste(refFiles, "/EPICArray/SNPProbes_McCartney.txt", sep = ""), stringsAsFactors = FALSE, header = TRUE)
crosshyb2<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM1.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes2<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM4.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes3<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM5.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes4<-read.csv(paste(refFiles, "/EPICArray/Pidsley_SM6.csv", sep = ""), stringsAsFactors = FALSE)
snpProbes<-snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01),]
snpProbes2<-snpProbes2[which(snpProbes2$AF > 0.01),]
snpProbes3<-snpProbes3[which(snpProbes3$AF > 0.01),]
snpProbes4<-snpProbes4[which(snpProbes4$AF > 0.01),]

dist<-cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist<-apply(dist, 1, min)
snpProbes4<-snpProbes4[which(dist <=10),]

remove<-intersect(rownames(celltypenormbeta),unique(c(crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE)))

celltypenormbeta<-celltypenormbeta[!rownames(celltypenormbeta) %in% remove,]


## calc cell type means including Total
cellMeans<-NULL
for(each in unique(pheno$Cell.type)){
	cellMeans<-cbind(cellMeans, apply(celltypenormbeta[,which(pheno$Cell.type == each)], 1, mean))
}
colnames(cellMeans)<-unique(pheno$Cell.type)

## calc cell type SDs including Total
cellSDs<-NULL
for(each in unique(pheno$Cell.type)){
	cellSDs<-cbind(cellSDs, apply(celltypenormbeta[,which(pheno$Cell.type == each)], 1, sd))
}
colnames(cellSDs)<-unique(pheno$Cell.type)


## filter to just cell fractions
pheno<-pheno[pheno$Cell.type %in% c("Double -ve","NeuN +ve","Sox10 +ve"),]
celltypenormbeta<-celltypenormbeta[,pheno$Basename]
pheno$Cell.type<-factor(pheno$Cell.type)

## anova for cell type differences at individual positions
ctOut<-apply(celltypenormbeta, 1, anovaModel, pheno$Cell.type, pheno$Age, pheno$Sex)
ctOut<-t(ctOut)
colnames(ctOut)<-c("NeuN:Coeff", "NeuN:P", "Sox10:Coeff", "Sox10:P", "ANOVA:P")

## filter out sites with missing location info
ctOut<-ctOut[!is.na(newfData$pos),]
celltypenormbeta<-celltypenormbeta[rownames(ctOut),]
newfData <- newfData[rownames(ctOut), ] 
ctOut<-cbind(ctOut, newfData[,c("chr","pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island", "Islands_Name")])
save(ctOut, file="Analysis/CellType/ANOVABrainCellTypes.rdata")

## how many loci?
dat<-data.frame("ids" = rownames(ctOut), "pval" = ctOut[,"NeuN:P"], "chr" = newfData$chr, "pos" = newfData$pos)
clump.neun<-clump(dat, thres1 = 9e-28, thres2 = 9e-28)

## how many genes?



## bumphunter to call regions

## filter samples with missing Age
celltypenormbeta<-celltypenormbeta[,!is.na(pheno$Age)]
pheno<-pheno[!is.na(pheno$Age),]
## sort
celltypenormbeta<-celltypenormbeta[order(newfData$chr, newfData$pos),]
newfData <- newfData[rownames(celltypenormbeta), ]
ctOut<-ctOut[rownames(newfData),]

cl<-clusterMaker(newfData$chr, newfData$pos, maxGap = 500)
#segs <- getSegments(ctOut[,1], cl, cutoff=0.05, assumeSorted = FALSE)
#tabs<-regionFinder(ctOut[,1],newfData$chr, newfData$pos, cl, cutoff = 0.05)
nCor<-detectCores()
registerDoParallel(cores = nCor)

designMat<-model.matrix(formula( ~ pheno$Cell.type + pheno$Age + pheno$Sex))
## need to run twice to compare Neun+ to all others and then Sox10+ to all others
tab.neun<-bumphunter(celltypenormbeta, designMat, newfData$chr, newfData$pos, cl, coef = 2, cutoff = 0.1, nullMethod="bootstrap")
tab.sox10<-bumphunter(celltypenormbeta, designMat, newfData$chr, newfData$pos, cl, coef = 3,cutoff = 0.1, nullMethod="bootstrap")
save(tab.neun, tab.sox10, file="Analysis/CellType/BumpHunterBrainCellTypes.rdata")

## as bootstrapping takes too long/never seems to finish even with a handful of iterations look at overlap of significant DMPs in these regions

dmps<-GRanges(seqnames = ctOut$chr, strand = "*", ranges = IRanges(start = ctOut$pos, end = ctOut$pos))
mcols(dmps)<-ctOut[,1:5]
neun.regions<-GRanges(seqnames = tab.neun$table$chr, strand = "*", ranges = IRanges(start = tab.neun$table$start, end = tab.neun$table$end))
mcols(neun.regions)<-tab.neun$table[,-c(1:3)]
sox10.regions<-GRanges(seqnames = tab.sox10$table$chr, strand = "*", ranges = IRanges(start = tab.sox10$table$start, end = tab.sox10$table$end))
mcols(sox10.regions)<-tab.sox10$table[,-c(1:3)]

## filter to regions with at least 2 sites
neun.regions<-neun.regions[neun.regions$L > 1]
sox10.regions<-sox10.regions[sox10.regions$L > 1]

## calculate proportion of probes in bump from cluster
summary(neun.regions$L/neun.regions$clusterL)
summary(sox10.regions$L/sox10.regions$clusterL)

## for each region count significant DMPs
interDmpsNeun<-findOverlaps(dmps, neun.regions)
nDMPs<-table(subjectHits(interDmpsNeun)[which(mcols(dmps)[queryHits(interDmpsNeun), "NeuN:P"] < 9e-8)])
neun.regions$nDMPs<-as.numeric(nDMPs[match(1:length(neun.regions), names(nDMPs))])
neun.regions$nDMPs[is.na(neun.regions$nDMPs)]<-0

interDmpsSox10<-findOverlaps(dmps, sox10.regions)
nDMPs<-table(subjectHits(interDmpsSox10)[which(mcols(dmps)[queryHits(interDmpsSox10), "Sox10:P"] < 9e-8)])
sox10.regions$nDMPs<-as.numeric(nDMPs[match(1:length(sox10.regions), names(nDMPs))])
sox10.regions$nDMPs[is.na(sox10.regions$nDMPs)]<-0

summary(neun.regions$nDMPs)
summary(neun.regions$nDMPs/neun.regions$L)
length(which(neun.regions$nDMPs/neun.regions$L == 1))
summary(sox10.regions$nDMPs)
summary(sox10.regions$nDMPs/sox10.regions$L)
length(which(sox10.regions$nDMPs/sox10.regions$L == 1))



## annotate with genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
neunRegionAnno <- annotatePeak(neun.regions, tssRegion=c(-3000, 3000), TxDb=txdb)
sox10RegionAnno <- annotatePeak(sox10.regions, tssRegion=c(-3000, 3000), TxDb=txdb)
save(neunRegionAnno, sox10RegionAnno, file = "Analysis/CellType/AnnotatedBumps.rdata")


## plot examples select those with the most DMPs
plotCols<-brewer.pal(4, "Set1")
axTrack <- GenomeAxisTrack()

pdf("Analysis/CellType/NeuNDMRs.pdf", width = 8, height = 8)
for(i in order(neun.regions$nDMPs, decreasing = TRUE)[1:10]){
	start <- tab.neun$table$start[i]
	end <- tab.neun$table$end[i]
	chr <- tab.neun$table$chr[i]

	## add in window around region
	windowSize<-end-start
	start<-start-windowSize
	end<-end+windowSize
	probesInCluster<-rownames(newfData)[which(newfData$chr == chr & newfData$pos <= end+1000000 & newfData$pos >= start-1000000)] ### add in extra sites so that plot starts and finishes outside of plotting region

	idxTrack <- IdeogramTrack(genome="hg19", chromosome=chr)
	refGenes <- UcscTrack(genome="hg19", chromosome=chr, table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=end, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
							symbol="name2", transcript="name", strand="strand", fill="#8282d2",
							  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")
	cpgIslands <- UcscTrack(genome="hg19", chromosome=chr, track="cpgIslandExt", from=start, to=end,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

	cellMeans.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), cellMeans[probesInCluster, ])
	dTrack1<-DataTrack(range = cellMeans.granges, genome = "hg19", type = "p", chromosome=chr, name = "DNAm Mean", col = plotCols, type = "smooth", groups = colnames(cellMeans), lwd = 2, baseline = seq(0,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	coeffs.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), ctOut[probesInCluster, c("NeuN:Coeff", "Sox10:Coeff")])
	pvals.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]),  -log10(ctOut[probesInCluster, c("NeuN:P")]))
	dTrack2<-DataTrack(range = coeffs.granges, genome = "hg19", type = "p", chromosome=chr, name = "Reg. Coeffs", col = plotCols[c(3,2)], type = "smooth", groups = c("NeuN", "Sox10"), lwd = 2, baseline = seq(-1,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	dTrack3<-DataTrack(range = pvals.granges, genome = "hg19", type = "p", chromosome=chr, name = "-logP", ylim = c(0, max(mcols(pvals.granges)[,1])), baseline = 0, col.baseline = "gray", lwd.baseline = 1)
	ht <- HighlightTrack(trackList=list(axTrack, refGenes,cpgIslands, dTrack1,dTrack2, dTrack3), start=tab.neun$table$start[i], end = tab.neun$table$end[i],chromosome=chr)

	## set plot to just inside
	plotTracks(list(idxTrack, ht), from=start, to=end, showTitle=TRUE)
}
dev.off()

pdf("Analysis/CellType/Sox10DMRs.pdf", width = 8, height = 8)
for(i in order(sox10.regions$nDMPs, decreasing = TRUE)[1:10]){
	start <- tab.sox10$table$start[i]
	end <- tab.sox10$table$end[i]
	chr <- tab.sox10$table$chr[i]

	## add in window around region
	windowSize<-end-start
	start<-start-windowSize
	end<-end+windowSize

	idxTrack <- IdeogramTrack(genome="hg19", chromosome=chr)
	refGenes <- UcscTrack(genome="hg19", chromosome=chr, table="ncbiRefSeq", track = 'NCBI RefSeq',from=start, to=end, trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
							symbol="name2", transcript="name", strand="strand", fill="#8282d2",
							  name="NCBI RefSeq", transcriptAnnotation="symbol", collapseTranscripts = "longest")
	cpgIslands <- UcscTrack(genome="hg19", chromosome=chr, track="cpgIslandExt", from=start, to=end,trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",shape="box", fill="#006400", name="CpG Islands")

	probesInCluster<-rownames(newfData)[which(newfData$chr == chr & newfData$pos <= end & newfData$pos >= start)]

	cellMeans.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), cellMeans[probesInCluster, ])
	dTrack1<-DataTrack(range = cellMeans.granges, genome = "hg19", type = "p", chromosome=chr, name = "DNAm Mean", col = plotCols, type = "smooth", groups = colnames(cellMeans), lwd = 2, baseline = seq(0,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	coeffs.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), ctOut[probesInCluster, c("NeuN:Coeff", "Sox10:Coeff")])
	pvals.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]),  -log10(ctOut[probesInCluster, c("Sox10:P")]))
	dTrack2<-DataTrack(range = coeffs.granges, genome = "hg19", type = "p", chromosome=chr, name = "Reg. Coeffs", col = plotCols[c(3,2)], type = "smooth", groups = c("NeuN", "Sox10"), lwd = 2, baseline = seq(-1,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	dTrack3<-DataTrack(range = pvals.granges, genome = "hg19", type = "p", chromosome=chr, name = "-logP", ylim = c(0, max(mcols(pvals.granges)[,1])), baseline = 0, col.baseline = "gray", lwd.baseline = 1)
	ht <- HighlightTrack(trackList=list(axTrack, refGenes,cpgIslands, dTrack1,dTrack2, dTrack3), start=tab.sox10$table$start[i], end = tab.sox10$table$end[i],chromosome=chr)

	## set plot to just inside
	plotTracks(list(idxTrack, ht), from=start, to=end, showTitle=TRUE)
}
dev.off()
