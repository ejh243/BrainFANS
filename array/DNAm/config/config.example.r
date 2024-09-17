## Study parameters
projectTitle<-"MRC Schizophrenia FANS samples"
processedBy<-"Complex Disease Epigenetic Group, University of Exeter Medical School"
tissueType<-"brain" # needs to be brain or blood to generate relevant cell composition estimates


## technical variables
techVar<-c("Sentrix_ID","Sentrix_Position")


## biological variables
bioVar<-c("Individual_ID","Cell_Type","Sex")


## QC thresholds
thresBS<-80
intenThres<-500
nvThres<-0.1
perMiss<-2


## Optional QC steps
sexCheck<-TRUE
snpCheck<-TRUE
ctCheck<-TRUE


## ctCheck variables
predDistinctCT<-c("NeuN+", "Sox10+", "IRF8+", "Triple-", "Total")
neunCT<-c("NeuN+", "SATB2+")
studentThres<-1.5
nSDThres<-3
