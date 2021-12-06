## this script is to be run prior to peak calling as assesses aligned reads
## this script combine multiple R packages and custom code to generate QC metrics to profile successful of ATAQ experiement
## in part based on this vignette: https://bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html
## the interpretation (provide in the rmd) is inspired by: https://galaxyproject.github.io/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#filtering-mapped-reads & https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3
## given the heavy computational burden, QC is split into batchs of 10 samples based on provided sample sheet
standardizeValues<-function(frag.len){
    x <- 1:1010
    frag.len <- frag.len[match(x, names(frag.len))]
    frag.len[is.na(frag.len)] <- 0
    y <- frag.len / sum(frag.len)
    y <- as.numeric(y)
    names(y) <- x
	return(y)
	}
	
calcNucleoProps<-function(fraglen.stand){
	freeSum<-sum(fraglen.stand[1:150])
	monoSum<-sum(fraglen.stand[151:300])
	diSum<-sum(fraglen.stand[301:450])
	triSum<-sum(fraglen.stand[451:600])
	otherSum<-sum(fraglen.stand[-c(1:600)])
	return(c(freeSum, monoSum, diSum, triSum, otherSum))
}

args <- commandArgs()
configFile<-args[6]
batchNum<-as.numeric(args[7]) ## nb starts from 0

source(configFile) ## contains paths to files

library(dplyr)
library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel") 
library(ATACseqQC)
library(Rsamtools)
library("phastCons100way.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38")
library(MotifDb)
library(ChIPpeakAnno)
register(DoparParam())
registered() 
bpparam("SerialParam")

setwd(dataDir) ## change to directory where aligned files (bam) and peaks (narrowPeaks) can be found ## will search for all within this folder

### Load Sample Sheet
pheno<-read.csv(sampleSheet, stringsAsFactors = FALSE)

## filter to required samples
index<-c(1:10)+(batchNum*10)
## if number of samples is not a function of ten adjust index
index<-index[index %in% 1:nrow(pheno)]

pheno<-pheno[index,]


## use ATACseqQC to calculate additional metrics
## some functions can be applied as lists
setwd(alignedDir)


## only makes sense to run on bam file from aligner (i.e no filtering!)
aFiles<-list.files(alignedDir, pattern = "_sorted.bam$", recursive = TRUE, full.names = TRUE)
aSampleNames<-gsub("_sorted.bam", "", basename(aFiles))

## subset out those in sample sheet
aFiles<-aFiles[match(paste(pheno$Project, pheno$Sample.Name, sep = "_"), aSampleNames)]
aSampleNames<-aSampleNames[match(pheno$Sample.Name, aSampleNames)]

## check all present
if(sum(is.na(aFiles))){
	print("Some samples not aligned")
	print("The following samples need aligning:")
	print(aFiles[is.na(aFiles)])
}

histDupReads<-lapply(aFiles,readsDupFreq) 
libComplexValues<-lapply(histDupReads, estimateLibComplexity)

## for this step find aligned QC'd reads chr1 only?
## only makes sense to run on bam file from aligner (i.e no filtering!)
aQCFiles<-list.files(alignedDir, pattern = "_depDup_q30.bam$", recursive = TRUE, full.names = TRUE)
aQCSampleNames<-gsub("_depDup_q30.bam", "", basename(aQCFiles))
aQCFiles<-aQCFiles[match(paste(pheno$Project, pheno$Sample.Name, sep = "_"), aQCSampleNames)]
aQCSampleNames<-aQCSampleNames[match(pheno$Sample.Name, aQCSampleNames)]

## remove NAs
aQCFiles<-aQCFiles[!is.na(aQCFiles)]
aQCSampleNames<-aQCSampleNames[!is.na(aQCSampleNames)]

fragSizeHist<-fragSizeDist(aQCFiles, aQCSampleNames)

## convert to ratios of nucleosomefree, mono, bi, tri etc


propNucleosomes<-lapply(lapply(fragSizeHist,standardizeValues),calcNucleoProps)
propNucleosomes<-matrix(data = unlist(propNucleosomes), ncol = 5, byrow = TRUE)

save(histDupReads, libComplexValues, fragSizeHist, propNucleosomes, file = paste0("QCoutput/ATACSeqQCObject_1_Batch", batchNum, ".rdata"))
#rm(list=c(histDupReads, libComplexValues, fragSizeHist, propNucleosomes))

## the following steps need to be run on each samples
#limit to chr1 to speed up processing
seqlev<-"chr1"


txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")


txs.chr1<- txs[seqnames(txs) %in% seqlev]
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
NTILE <- 51
dws <- ups <- 1010


## create GRanges object of genomic regions (i.e. chr1) to load
which <- GRanges(seqlev, IRanges(0, seqlengths(txs.chr1)[seqlev]))
#p1 <- ScanBamParam(what=scanBamWhat(), which=which)

## other need to be applied in loop below
pt.collate<-list(length = length(aQCFiles))
nfr.collate<-list(length = length(aQCFiles))
tsse.collate<-list(length = length(aQCFiles))
sigs.collate<-list(length = length(aQCFiles))
names(pt.collate)<-aQCSampleNames
names(nfr.collate)<-aQCSampleNames
names(tsse.collate)<-aQCSampleNames
names(sigs.collate)<-aQCSampleNames
nFragtype.collate<-matrix(data = NA, nrow = length(aQCFiles), ncol = 8)
rownames(nFragtype.collate)<-aQCSampleNames
for(i in 1:length(aQCFiles)){
	print(paste("Running Sample:", aQCSampleNames[i]))
	bamfile<-aQCFiles[i]
	#bamfile<-gsub("\\.bam", "_chr1\\.bam", bamfile)
	bamDat <- readBamFile(bamfile, asMates=TRUE, bigFile=TRUE, which = which)
	
	## adjust the start sites of the reads, by default, all reads aligning to the positive strand are offset by +4bp, and all reads aligning to the negative strand are offset by -5bp.
	## save output for peak calling
	bamDat.shift <- shiftGAlignmentsList(bamDat)
	## save shifted reads as useful for downstream analysis
	#shiftedBamfile <- gsub("\\.bam","_shifted\\.bam",  bamfile)
	#export(bamDat.shift, shiftedBamfile)
	## calculate promotor/transcript score
	pt.collate[[i]] <- PTscore(bamDat.shift, txs)
	## calculate Nucleosome Free Regions (NFR) score
	nfr.collate[[i]] <- NFRscore(bamDat.shift, txs)
	# Transcription Start Site (TSS) Enrichment Score
	tsse.collate[[i]] <- TSSEscore(bamDat.shift, txs)

	# Split reads into nucleosome free,mono, dinucleo, etc ...
	## NB time consuming step, quicker if exclude conservation matching
	#objs <- splitGAlignmentsByCut(bamDat.shift, outPath=paste0("splitted", batchNum),
    #             txs=txs.chr1, genome=genome,
     #
	objs <- splitGAlignmentsByCut(bamDat.shift, outPath=paste0("splitted", batchNum),
                 txs=txs.chr1, genome=genome)
	
	# count how many assigned to each class
	nFragtype.collate[i,]<-unlist(lapply(objs, length))
	
	bamfiles <- file.path(paste0("splitted", batchNum),
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))
					 
	## Plot the cumulative percentage of tag allocation in nucleosome-free 
	## and mononucleosome bam files.
	#cumulativePercentage(bamfiles[1:2], as(seqinfo(Hsapiens)[seqlev], "GRanges"))
	
	## heatmap and coverage curve for nucleosome positions
	## estmate library size for normalisation
	librarySize <- estLibSize(bamfiles[1:3])
	
	## calculate the signals around TSSs.
	sigs.collate[[i]] <- enrichedFragments(gal=objs[c("NucleosomeFree", 
										 "mononucleosome",
										 "dinucleosome",
										 "trinucleosome")], 
							  TSS=TSS,
							  librarySize=librarySize,
							  seqlev=seqlev,
							  TSS.filter=0.5,
							  n.tile = NTILE,
							  upstream = ups,
							  downstream = dws)
						  
}


save(pt.collate, nfr.collate, sigs.collate, tsse.collate,nFragtype.collate,  file = paste0("QCoutput/ATACSeqQCObject_2_Batch", batchNum, ".rdata"))

