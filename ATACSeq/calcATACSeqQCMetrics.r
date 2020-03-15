## this script is to be run prior to peak calling as assesses aligned reads
## this script combine multiple R packages and custom code to generate QC metrics to profile successful of ATAQ experiement
## in part based on this vignette: https://bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html
## the interpretation (provide in the rmd) is inspired by: https://galaxyproject.github.io/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#filtering-mapped-reads & https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3

source("config.r") ## contains paths to files

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

### Create Sample Sheet
bamReads<-list.files(".", pattern = "_depDup_q30.bam", recursive = TRUE)
bamReads<-bamReads[grep("bai", bamReads, invert = TRUE)]
bamIDs<-unlist(lapply(strsplit(bamReads, "/"), tail, n = 1))
bamIDs<-gsub("_trimmed_depDup_q30.bam", "", bamIDs)

tissue<-unlist(lapply(strsplit(bamIDs, "_"),tail, n = 1))

pe<-"Paired"

sampleSheet<-data.frame(SampleID = bamIDs, Tissue=tissue, Factor="ATAC", ReadType = pe, bamReads = bamReads, stringsAsFactors = FALSE)

#write.csv(sampleSheet, "SampleSheetForATACSeqQC.csv")

## load uniqueness info
histFiles<-list.files(pattern = "_hist.txt", recursive = TRUE)
fileNames<-gsub("_trimmed_r1.fq.gz_hist.txt", "", histFiles)
fileNames<-gsub("_trimmed_hist.txt", "", fileNames)
fileNames<-gsub("11_trimmed/", "", fileNames)

hist.data<-read.table(histFiles[1])
if(ncol(hist.data) == 7){
	colnames(hist.data)<-c("count","first","rand","first_cnt","pair_cnt","avg_quality","perfect_prob")
} else{
	if(ncol(hist.data) == 17){
	colnames(hist.data)<-c("count","first","rand","r1_first","r1_rand","r2_first","r2_rand","pair","first_cnt","rand_cnt","r1_first_cnt","r1_rand_cnt","r2_first_cnt","r2_rand_cnt","pair_cnt","avg_quality","perfect_prob")
	
	}
}

hist.data<-hist.data[,c("count", "pair_cnt")]
for(each in histFiles[-1]){

	tmp.dat<-read.table(each)
	if(ncol(tmp.dat) == 7){
		colnames(tmp.dat)<-c("count","first","rand","first_cnt","pair_cnt","avg_quality","perfect_prob")
		tmp.dat<-tmp.dat[,c("count", "pair_cnt")]
	} else{
		if(ncol(tmp.dat) == 17){
		colnames(tmp.dat)<-c("count","first","rand","r1_first","r1_rand","r2_first","r2_rand","pair","first_cnt","rand_cnt","r1_first_cnt","r1_rand_cnt","r2_first_cnt","r2_rand_cnt","pair_cnt","avg_quality","perfect_prob")
		## take count unique as mean of r1 and r2_first
		tmp.dat$pair_cnt<-apply(tmp.dat[, c("r1_first_cnt", "r2_first_cnt")], 1, mean)
		tmp.dat<-tmp.dat[,c("count", "pair_cnt")]
	}
	}

	
	hist.data<-full_join(hist.data, tmp.dat, by = "count")
}


colnames(hist.data)[-1]<-fileNames

## use ATACseqQC to calculate additional metrics


## some functions can be applied as lists
histDupReads<-lapply(gsub("_depDup_q30", "_sorted", sampleSheet$bamReads),readsDupFreq) ## only makes sens to run on bam file from aligner (i.e no filtering!)
libComplexValues<-lapply(histDupReads, estimateLibComplexity)
fragSizeHist<-fragSizeDist(sampleSheet$bamReads, sampleSheet$SampleID)

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
seqlev<-"chr1"

txs.chr1<- txs[seqnames(txs) %in% seqlev]
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
NTILE <- 51
dws <- ups <- 1010

## other need to be applied in loop below
pt.collate<-list(length = nrow(sampleSheet))
nfr.collate<-list(length = nrow(sampleSheet))
tsse.collate<-list(length = nrow(sampleSheet))
sigs.collate<-list(length = nrow(sampleSheet))

for(i in 1:nrow(sampleSheet)){
	bamfile<-sampleSheet$bamReads[i]
	bamDat <- readBamFile(bamfile, asMates=TRUE, bigFile=TRUE)
	shiftedBamfile <- gsub("alignedData","shifted",  bamfile)
	
	## adjust the start sites of the reads, by default, all reads aligning to the positive strand are offset by +4bp, and all reads aligning to the negative strand are offset by -5bp.
	## save output for peak calling
	bamDat.shift <- shiftGAlignmentsList(bamDat)
	## calculate promotor/transcript score
	pt.collate[[i]] <- PTscore(bamDat.shift, txs)
	## calculate Nucleosome Free Regions (NFR) score
	nfr.collate[[i]] <- NFRscore(bamDat.shift, txs)
	# Transcription Start Site (TSS) Enrichment Score
	tsse.collate[[i]] <- TSSEscore(bamDat.shift, txs)
	# Split reads NB time consuming step so filter to single chr? ## also quick if exclude conservation matching
	objs <- splitGAlignmentsByCut(bamDat.shift, outPath="splitted",
                 txs=txs.chr1, genome=genome,
                 conservation=phastCons100way.UCSC.hg38)
				 
	bamfiles <- file.path("splitted",
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))
	## Plot the cumulative percentage of tag allocation in nucleosome-free 
	## and mononucleosome bam files.
	#cumulativePercentage(bamfiles[1:2], as(seqinfo(Hsapiens)["chr1"], "GRanges")) # couldn't get this to work
	
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


save(hist.data, histDupReads, libComplexValues, fragSizeHist, pt.collate, nfr.collate, sigs.collate, tsse.collate, file = "ATACSeqQCObject_1.rdata")

