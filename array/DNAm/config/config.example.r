## Study parameters
projectTitle<-"MRC Schizophrenia FANS samples"
processedBy<-"Complex Disease Epigenetic Group, University of Exeter Medical School"
arrayVersion<-"Illumina EPIC microarray"
tissueType<-"brain" # needs to be brain or blood to generate relevant cell composition estimates
cellSorted<-"FALSE"

## technical variables

techVar<-c("Sentrix_ID","Sentrix_Position")

## biological variables

bioVar<-c("Individual_ID","Cell_Type","Sex")

predDistinctCT<-c("NeuN+", "Sox10+", "IRF8+", "Triple-", "Total")
neunCT<-c("NeuN+", "SATB2+")
otherCT<-c("Double-", "IRF8+", "SATB2", "Sox10+", "Triple-")



## QC thresholds
thresBS<-80
intenThres<-500
nvThres<-0.1
sexCheck<-TRUE
snpCheck<-TRUE
ctCheck<-TRUE
perMiss<-2
studentThres<-1.5
nSDThres<-3
pnthres<-0.05
perc<-1
