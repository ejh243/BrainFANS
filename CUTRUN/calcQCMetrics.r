#source("config.r") ## contains paths to files
## need to extend to run additional peak sets called by CUT&RUNtools pipeline
## use to compare merits of different post-processed files and peak sets. 


alignedDir<-"AlignedData/"
peakDir<-"CalledPeaks/"


library(ChIPQC)
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library("BiocParallel") 
library(reshape)


### Create Sample Sheet
bamReads<-list.files(paste0(alignedDir, "sorted/"), pattern = ".bam$", recursive = TRUE)
bamIDs<-gsub(".bam", "", bamReads)
sampleInfo<-matrix(data = unlist(strsplit(bamIDs, "_")), ncol = 4, byrow = TRUE)
colnames(sampleInfo)<-c("SequencingProject", "Mark", "Individual.ID", "CellType")
sampleInfo<-as.data.frame(sampleInfo)

## use narrow peak calls 
narrowPeaks<-list.files(paste0(peakDir, "macs2.narrow.aug18"), pattern = ".narrowPeak", recursive = TRUE)
peakIDs<-gsub("_peaks.narrowPeak", "", narrowPeaks)

sampleSheet<-data.frame(SampleID = bamIDs, Tissue=sampleInfo$CellType, Factor=sampleInfo$Mark, Replicate=1, 
ReadType = "Paired", bamReads = paste(alignedDir, "sorted", bamReads, sep = "/"), 
Peaks = paste(peakDir, "macs2.narrow.aug18", narrowPeaks[match(bamIDs, peakIDs)], sep = "/"), stringsAsFactors = FALSE)

sampleSheet$Peaks[endsWith(sampleSheet$Peaks, "NA")]<-NA

write.csv(sampleSheet, "QCmetrics/SampleSheetForCUTRUN.csv") ## overwrites existing samples sheet

dat<-ChIPQC(sampleSheet, consensus = FALSE, chromosomes = NULL, annotation = "hg38")

save(dat, file = "QCmetrics/ChIPQCMetrics.rdata")


par(mfrow = c(2,2))
par(mar = c(5,5,1,1))
#plotPeakProfile(dat, facetBy = FALSE) + theme(text = element_text(size=20))
PeakSignal <- averagepeaksignal(dat)
  if(all(!is.na(PeakSignal))){
    Width <- seq(-nrow(PeakSignal)/2,nrow(PeakSignal)/2)[-(nrow(PeakSignal)/2+1)]
    Window <- length(Width)
    PSDataFrame <- data.frame(Width,PeakSignal)
    PSDataFrame <- melt(PSDataFrame,id.vars=c("Width"))
    colnames(PSDataFrame)<-c("Distance", "Sample", "Signal")
	PSDataFrame$Mark<-sampleSheet$Factor[match(PSDataFrame$Sample, paste0("X", sampleSheet$SampleID))]
	PSDataFrame$CellType<-sampleSheet$Tissue[match(PSDataFrame$Sample, paste0("X", sampleSheet$SampleID))]

    Plot <- ggplot(PSDataFrame,aes(x=Distance,y=Signal, group=Sample))+
      geom_line(size=1.3, aes(col = Mark))+xlim(-Window/2,Window/2)+ylab("Signal")+
      theme(axis.title.y=element_text(angle=0)) + facet_wrap(~CellType)
    Plot
  }else{
    stop("No average signal available for sample")
  }
  
 plotRap(dat, facet = FALSE) + theme(text = element_text(size=20), axis.title = element_text(size = 10, angle = 25))
 
 par(mar = c(5,5,1,1))
y_lim<-c(0,max(frip(dat)*100, na.rm = TRUE)*1.1)
barplot(frip(dat)*100, ylim = y_lim, ylab = "% Reads in Peaks", names.arg = sampleSheet$IID, 
        cex.names = 2, cex.axis = 2, cex.lab = 2)

