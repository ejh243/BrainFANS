## Study parameters
projectTitle<-"MRC Schizophrenia FANS samples"
processedBy<-"Complex Disease Epigenetic Group, University of Exeter Medical School"
tissueType<-"brain" # needs to be brain or blood to generate relevant cell composition estimates
arrayType<-"V2" # must be one of "450K" / "V1" / "V2"


## project variables
projVar<-c("Cell_Type", "Sex") # Used to colour plots


## QC thresholds
thresBS<-80
intenThres<-500
nvThres<-0.1
perMiss<-2


## Optional QC steps
sexCheck<-TRUE
snpCheck<-TRUE
ctCheck<-TRUE
extractScanDate<-FALSE # Adds scan date to gds files


## ctCheck variables
predDistinctCT<-c("NeuN+", "Sox10+", "IRF8+", "Triple-", "Total")
neunCT<-c("NeuN+", "SATB2+")
studentThres<-1.5
nSDThres<-3
