anovaModel<-function(row, celltype, age,sex){
	model<-lm(row ~ celltype + age + sex)
	null<-lm(row ~ age + sex)
	tmp<-summary(model)$coefficients
	regCoeffs<-c(t(tmp[c("celltypeNeuN +ve", "celltypeSox10 +ve", "celltypeIRF8 +ve"), c(1,4)]), anova(model,null)[2,6])
	return(regCoeffs)
}

library(ChIPseeker)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(doParallel)


### test cell comp differences
source("rmdConfig.mrc")
setwd(dataDir)
load(normData)

## remove rs probes
celltypenormbeta<-celltypenormbeta[grep("rs", rownames(celltypenormbeta), invert = TRUE),]


## calc cell type means including Total
cellMeans<-NULL
for(each in sort(unique(pheno$Cell.type))){
	cellMeans<-cbind(cellMeans, apply(celltypenormbeta[,which(pheno$Cell.type == each)], 1, mean))
}
colnames(cellMeans)<-sort(unique(pheno$Cell.type))

## calc cell type SDs including Total
cellSDs<-NULL
for(each in sort(unique(pheno$Cell.type))){
	cellSDs<-cbind(cellSDs, apply(celltypenormbeta[,which(pheno$Cell.type == each)], 1, sd))
}
colnames(cellSDs)<-sort(unique(pheno$Cell.type))


## filter to just cell fractions
pheno<-pheno[pheno$Cell.type %in% c("Double -ve","NeuN +ve","Sox10 +ve", "IRF8 +ve"),]
celltypenormbeta<-celltypenormbeta[,pheno$Basename]
pheno$Cell.type<-factor(pheno$Cell.type)

## anova for cell type differences at individual positions
ctOut<-apply(celltypenormbeta, 1, anovaModel, pheno$Cell.type, pheno$Age, pheno$Sex)
ctOut<-t(ctOut)
colnames(ctOut)<-c("NeuN:Coeff", "NeuN:P", "Sox10:Coeff", "Sox10:P",  "IRF8:Coeff", "IRF8:P", "ANOVA:P")

annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))

newfData<-newfData[rownames(celltypenormbeta),]

## filter out sites with missing location info
ctOut<-ctOut[!is.na(newfData$pos),]
celltypenormbeta<-celltypenormbeta[rownames(ctOut),]
newfData <- newfData[rownames(ctOut), ] 
ctOut<-cbind(ctOut, newfData[,c("chr","pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island", "Islands_Name")])
save(ctOut, file="Analysis/CellType/ANOVABrainCellTypes.rdata")

length(which(ctOut[,7] < 9e-8))
length(which(ctOut[,2] < 9e-8))
length(which(ctOut[,4] < 9e-8))
length(which(ctOut[,6] < 9e-8))

pdf("Analysis/CellType/QQplot.pdf")
par(mfrow = c(2,2))
qq(ctOut[,7], main = "Any")
qq(ctOut[,2], main = "NeuNvsDoubleNeg")
qq(ctOut[,4], main = "Sox10vsDoubleNeg")
qq(ctOut[,6], main = "IRF8vsDoubleNeg")
dev.off()

index<-which(ctOut[,7] < 9e-8)

## look at direction of effect
pdf("Analysis/CellType/HistDMPsMeanDifferences.pdf", width = 12, height = 4)
par(mfrow = c(1,3))
hist(ctOut[which(ctOut[,2] < 9e-8),1], xlab = "Mean Difference", main = "NeuNvsDoubleNegative", ylab = "nDMPs")
signTab<-table(sign(ctOut[which(ctOut[,2] < 9e-8),1]))
signTest<-binom.test(signTab[1], sum(signTab))
text(x =  par("usr")[2],y =  par("usr")[4], pos = 2, offset = -1, paste(signif(signTab[2]/sum(signTab)*100,3), "% DMPs hypermethylated.\nP =", signif(signTest$p.value, 3)), xpd= TRUE)
hist(ctOut[which(ctOut[,4] < 9e-8),3], xlab = "Mean Difference", main = "Sox10vsDoubleNegative", ylab = "nDMPs")
signTab<-table(sign(ctOut[which(ctOut[,4] < 9e-8),3]))
signTest<-binom.test(signTab[1], sum(signTab))
text(x =  par("usr")[2],y =  par("usr")[4], pos = 2, offset = -1, paste(signif(signTab[2]/sum(signTab)*100,3), "% DMPs hypermethylated.\nP =", signif(signTest$p.value, 3)), xpd= TRUE)
hist(ctOut[which(ctOut[,6] < 9e-8),5], xlab = "Mean Difference", main = "IRF8vsDoubleNegative", ylab = "nDMPs")
signTab<-table(sign(ctOut[which(ctOut[,6] < 9e-8),5]))
signTest<-binom.test(signTab[1], sum(signTab))
text(x =  par("usr")[2],y =  par("usr")[4], pos = 2, offset = -1, paste(signif(signTab[2]/sum(signTab)*100,3), "% DMPs hypermethylated.\nP =", signif(signTest$p.value, 3)), xpd= TRUE)
dev.off()

library(VennDiagram)
library(minfi)
  
# Chart
venn.diagram(
  x = list(rownames(ctOut)[which(ctOut[,2] < 9e-8)], rownames(ctOut)[which(ctOut[,4] < 9e-8)], rownames(ctOut)[which(ctOut[,6] < 9e-8)]),
  category.names = c("NeuNvsDoubleNeg" , "Sox10vsDoubleNeg" , "IRF8vsDoubleNeg"),
  filename = 'Analysis/CellType/VennDiagramCelltypeDMPs.png',
  output=TRUE
)

## for probes that are DMPs in two cell types, same direction and magnitude?

Neun.ES<-cut(ctOut[,1], seq(-1,1,0.05))
Sox10.ES<-cut(ctOut[,3], seq(-1,1,0.05))
Irf8.ES<-cut(ctOut[,5], seq(-1,1,0.05))

pdf("Analysis/CellType/HeatmapSharedDMPsBetweenCellTypes.pdf")
index<-which(ctOut[,2] < 9e-8 & ctOut[,4] < 9e-8)
heatmap(table(Neun.ES[index], Sox10.ES[index]), Rowv = NA, Colv = NA, xlab = "NeuN+ve vs Double-ve", ylab = "Sox10+ve vs Double-ve")
index<-which(ctOut[,2] < 9e-8 & ctOut[,6] < 9e-8)
heatmap(table(Neun.ES[index], Irf8.ES[index]), Rowv = NA, Colv = NA, xlab = "NeuN+ve vs Double-ve", ylab = "Irf8+ve vs Double-ve")
index<-which(ctOut[,4] < 9e-8 & ctOut[,6] < 9e-8)
heatmap(table(Sox10.ES[index], Irf8.ES[index]), Rowv = NA, Colv = NA, xlab = "Sox10+ve vs Double-ve", ylab = "Irf8+ve vs Double-ve")
dev.off()



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
tab.neun<-bumphunter(celltypenormbeta, designMat, newfData$chr, newfData$pos, cl, coef = 3, cutoff = 0.1, nullMethod="bootstrap")
tab.irf8<-bumphunter(celltypenormbeta, designMat, newfData$chr, newfData$pos, cl, coef = 2,cutoff = 0.1, nullMethod="bootstrap")
tab.sox10<-bumphunter(celltypenormbeta, designMat, newfData$chr, newfData$pos, cl, coef = 4,cutoff = 0.1, nullMethod="bootstrap")
save(tab.neun, tab.sox10, tab.irf8,file="Analysis/CellType/BumpHunterBrainCellTypes.rdata")

## as bootstrapping takes too long/never seems to finish even with a handful of iterations look at overlap of significant DMPs in these regions

dmps<-GRanges(seqnames = ctOut$chr, strand = "*", ranges = IRanges(start = ctOut$pos, end = ctOut$pos))
mcols(dmps)<-ctOut[,1:7]
neun.regions<-GRanges(seqnames = tab.neun$table$chr, strand = "*", ranges = IRanges(start = tab.neun$table$start, end = tab.neun$table$end))
mcols(neun.regions)<-tab.neun$table[,-c(1:3)]
sox10.regions<-GRanges(seqnames = tab.sox10$table$chr, strand = "*", ranges = IRanges(start = tab.sox10$table$start, end = tab.sox10$table$end))
mcols(sox10.regions)<-tab.sox10$table[,-c(1:3)]
irf8.regions<-GRanges(seqnames = tab.irf8$table$chr, strand = "*", ranges = IRanges(start = tab.irf8$table$start, end = tab.irf8$table$end))
mcols(irf8.regions)<-tab.irf8$table[,-c(1:3)]

## filter to regions with at least 2 sites
neun.regions<-neun.regions[neun.regions$L > 1]
sox10.regions<-sox10.regions[sox10.regions$L > 1]
irf8.regions<-irf8.regions[irf8.regions$L > 1]

## calculate proportion of probes in bump from cluster
summary(neun.regions$L/neun.regions$clusterL)
summary(sox10.regions$L/sox10.regions$clusterL)
summary(irf8.regions$L/irf8.regions$clusterL)

## for each region count significant DMPs
interDmpsNeun<-findOverlaps(dmps, neun.regions)
nDMPs<-table(subjectHits(interDmpsNeun)[which(mcols(dmps)[queryHits(interDmpsNeun), "NeuN:P"] < 9e-8)])
neun.regions$nDMPs<-as.numeric(nDMPs[match(1:length(neun.regions), names(nDMPs))])
neun.regions$nDMPs[is.na(neun.regions$nDMPs)]<-0

interDmpsSox10<-findOverlaps(dmps, sox10.regions)
nDMPs<-table(subjectHits(interDmpsSox10)[which(mcols(dmps)[queryHits(interDmpsSox10), "Sox10:P"] < 9e-8)])
sox10.regions$nDMPs<-as.numeric(nDMPs[match(1:length(sox10.regions), names(nDMPs))])
sox10.regions$nDMPs[is.na(sox10.regions$nDMPs)]<-0

interDmpsIrf8<-findOverlaps(dmps, irf8.regions)
nDMPs<-table(subjectHits(interDmpsIrf8)[which(mcols(dmps)[queryHits(interDmpsIrf8), "IRF8:P"] < 9e-8)])
irf8.regions$nDMPs<-as.numeric(nDMPs[match(1:length(irf8.regions), names(nDMPs))])
irf8.regions$nDMPs[is.na(irf8.regions$nDMPs)]<-0


summary(neun.regions$nDMPs)
summary(neun.regions$nDMPs/neun.regions$L)
length(which(neun.regions$nDMPs/neun.regions$L == 1))
summary(sox10.regions$nDMPs)
summary(sox10.regions$nDMPs/sox10.regions$L)
length(which(sox10.regions$nDMPs/sox10.regions$L == 1))
summary(irf8.regions$nDMPs)
summary(irf8.regions$nDMPs/irf8.regions$L)
length(which(irf8.regions$nDMPs/irf8.regions$L == 1))

## summarise size of bumps
summary(width(neun.regions))
summary(width(sox10.regions))
summary(width(irf8.regions))
summary(width(neun.regions[which(neun.regions$nDMPs > 0),]))
summary(width(sox10.regions[which(sox10.regions$nDMPs > 0),]))
summary(width(irf8.regions[which(irf8.regions$nDMPs > 0),]))
summary(width(neun.regions[which(neun.regions$nDMPs/neun.regions$L == 1),]))
summary(width(sox10.regions[which(sox10.regions$nDMPs/sox10.regions$L == 1),]))
summary(width(irf8.regions[which(irf8.regions$nDMPs/irf8.regions$L == 1),]))

## identify regions specific to each cell type
neun.specific<-neun.regions[unique(subjectHits(findOverlaps(setdiff(neun.regions, reduce(c(sox10.regions, irf8.regions))), neun.regions))),]
sox10.specific<-sox10.regions[unique(subjectHits(findOverlaps(setdiff(sox10.regions,reduce(c(neun.regions, irf8.regions))), sox10.regions))),]
irf8.specific<-irf8.regions[unique(subjectHits(findOverlaps(setdiff(irf8.regions,reduce(c(neun.regions, sox10.regions))), irf8.regions))),]

## noticed some specific regions are consecuative
neun.conseq<-unique(c(neun.specific[queryHits(findOverlaps(neun.specific, sox10.specific)),], neun.specific[queryHits(findOverlaps(neun.specific, irf8.specific)),]))
sox10.conseq<-unique(c(sox10.specific[subjectHits(findOverlaps(neun.specific, sox10.specific)),],sox10.specific[subjectHits(findOverlaps(irf8.specific, sox10.specific)),]))
irf8.conseq<-unique(c(irf8.specific[subjectHits(findOverlaps(neun.specific, irf8.specific)),],irf8.specific[subjectHits(findOverlaps(sox10.specific, irf8.specific)),]))


## how many DMPs are not located within 

singleProbeCl<-names(table(cl))[which(table(cl) == 1)]
table(unique(cl[which(ctOut[,2] < 9e-8)]) %in% singleProbeCl)
table(unique(cl[which(ctOut[,4] < 9e-8)]) %in% singleProbeCl)
table(unique(cl[which(ctOut[,6] < 9e-8)]) %in% singleProbeCl)


## annotate with genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
neunRegionAnno <- annotatePeak(neun.regions, tssRegion=c(-3000, 3000), TxDb=txdb)
sox10RegionAnno <- annotatePeak(sox10.regions, tssRegion=c(-3000, 3000), TxDb=txdb)
irf8RegionAnno <- annotatePeak(irf8.regions, tssRegion=c(-3000, 3000), TxDb=txdb)
save(neunRegionAnno, sox10RegionAnno, irf8RegionAnno, file = "Analysis/CellType/AnnotatedBumps.rdata")


## plot examples select those with the most DMPs
plotCols<-c("darkgreen", "darkblue", "darkmagenta", "deeppink", "darkgray") ## assumes celltypes are ordered alphabetically

axTrack <- GenomeAxisTrack()

pdf("Analysis/CellType/NeuNSpecificDMRs.pdf", width = 8, height = 8)
for(i in order(neun.specific$nDMPs, decreasing = TRUE)[1:10]){
	start <- as.numeric(start(neun.specific)[i])
	end <- as.numeric(end(neun.specific)[i])
	chr <- as.character(seqnames(neun.specific)[i])
	
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
	dTrack1<-DataTrack(range = cellMeans.granges, genome = "hg19", type = "p", chromosome=chr, name = "DNAm Mean", col = plotCols[1:4], type = "smooth", groups = c("Double -ve", "IRF8+ve", "NeuN+ve", "Sox10+ve"), lwd = 2, baseline = seq(0,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	coeffs.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), ctOut[probesInCluster, c("IRF8:Coeff", "NeuN:Coeff", "Sox10:Coeff")])
	pvals.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), "IRF8:log10P" = -log10(ctOut[probesInCluster, c("IRF8:P")]),"NeuN:log10P" =  -log10(ctOut[probesInCluster, c("NeuN:P")]), "Sox10:log10P" = -log10(ctOut[probesInCluster, c("Sox10:P")]))
	dTrack2<-DataTrack(range = coeffs.granges, genome = "hg19", type = "p", chromosome=chr, name = "Reg. Coeffs", col = plotCols[c(2:4)], type = "smooth", groups = c("IRF8+ve", "NeuN+ve", "Sox10+ve"), lwd = 2, baseline = seq(-1,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	dTrack3<-DataTrack(range = pvals.granges, genome = "hg19", type = "p", chromosome=chr, name = "-logP", ylim = c(0, max(mcols(pvals.granges)[,1])), groups = c("IRF8+ve", "NeuN+ve", "Sox10+ve"), baseline = 0, col.baseline = "gray", lwd.baseline = 1, col = plotCols[c(2:4)])
	ht <- HighlightTrack(trackList=list(axTrack, refGenes,cpgIslands, dTrack1,dTrack2, dTrack3), start=as.numeric(start(neun.specific)[i]), end = as.numeric(end(neun.specific)[i]),chromosome=chr)

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
	dTrack1<-DataTrack(range = cellMeans.granges, genome = "hg19", type = "p", chromosome=chr, name = "DNAm Mean", col = plotCols[1:4], type = "smooth", groups = c("Double -ve", "IRF8+ve", "NeuN+ve", "Sox10+ve"), lwd = 2, baseline = seq(0,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	coeffs.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), ctOut[probesInCluster, c("IRF8:Coeff", "NeuN:Coeff", "Sox10:Coeff")])
	pvals.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), "IRF8:log10P" = -log10(ctOut[probesInCluster, c("IRF8:P")]),"NeuN:log10P" =  -log10(ctOut[probesInCluster, c("NeuN:P")]), "Sox10:log10P" = -log10(ctOut[probesInCluster, c("Sox10:P")]))
	dTrack2<-DataTrack(range = coeffs.granges, genome = "hg19", type = "p", chromosome=chr, name = "Reg. Coeffs", col = plotCols[c(2:4)], type = "smooth", groups = c("IRF8+ve", "NeuN+ve", "Sox10+ve"), lwd = 2, baseline = seq(-1,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	dTrack3<-DataTrack(range = pvals.granges, genome = "hg19", type = "p", chromosome=chr, name = "-logP", ylim = c(0, max(mcols(pvals.granges)[,1])), groups = c("IRF8+ve", "NeuN+ve", "Sox10+ve"), baseline = 0, col.baseline = "gray", lwd.baseline = 1, col = plotCols[c(2:4)])
	ht <- HighlightTrack(trackList=list(axTrack, refGenes,cpgIslands, dTrack1,dTrack2, dTrack3), start=tab.sox10$table$start[i], end = tab.sox10$table$end[i],chromosome=chr)

	## set plot to just inside
	plotTracks(list(idxTrack, ht), from=start, to=end, showTitle=TRUE)
}
dev.off()


pdf("Analysis/CellType/IRF8DMRs.pdf", width = 8, height = 8)
for(i in order(irf8.regions$nDMPs, decreasing = TRUE)[1:10]){
	start <- tab.irf8$table$start[i]
	end <- tab.irf8$table$end[i]
	chr <- tab.irf8$table$chr[i]

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
	dTrack1<-DataTrack(range = cellMeans.granges, genome = "hg19", type = "p", chromosome=chr, name = "DNAm Mean", col = plotCols[1:4], type = "smooth", groups = c("Double -ve", "IRF8+ve", "NeuN+ve", "Sox10+ve"), lwd = 2, baseline = seq(0,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	coeffs.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), ctOut[probesInCluster, c("IRF8:Coeff", "NeuN:Coeff", "Sox10:Coeff")])
	pvals.granges<-GRanges(seqnames = newfData[probesInCluster, "chr"], strand = "*", ranges = IRanges(start = newfData[probesInCluster, "pos"], end = newfData[probesInCluster, "pos"]), "IRF8:log10P" = -log10(ctOut[probesInCluster, c("IRF8:P")]),"NeuN:log10P" =  -log10(ctOut[probesInCluster, c("NeuN:P")]), "Sox10:log10P" = -log10(ctOut[probesInCluster, c("Sox10:P")]))
	dTrack2<-DataTrack(range = coeffs.granges, genome = "hg19", type = "p", chromosome=chr, name = "Reg. Coeffs", col = plotCols[c(2:4)], type = "smooth", groups = c("IRF8+ve", "NeuN+ve", "Sox10+ve"), lwd = 2, baseline = seq(-1,1,0.2), col.baseline = "gray", lwd.baseline = 1)
	dTrack3<-DataTrack(range = pvals.granges, genome = "hg19", type = "p", chromosome=chr, name = "-logP", ylim = c(0, max(mcols(pvals.granges)[,1])), groups = c("IRF8+ve", "NeuN+ve", "Sox10+ve"), baseline = 0, col.baseline = "gray", lwd.baseline = 1, col = plotCols[c(2:4)])
	ht <- HighlightTrack(trackList=list(axTrack, refGenes,cpgIslands, dTrack1,dTrack2, dTrack3), start=tab.irf8$table$start[i], end = tab.irf8$table$end[i],chromosome=chr)

	## set plot to just inside
	plotTracks(list(idxTrack, ht), from=start, to=end, showTitle=TRUE)
}
dev.off()

