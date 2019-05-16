## Written by Eilis
## Design 5hmC array plate layout profiling
## 1 ID per chip 
## randomise the order of sample type but keep OXBS and BS next to each other

setwd("U:/Projects/Schizophrenia/Brain/PlateDesigns")
sampleSheet<-read.table("BSOxBSFACSsort_Chip loc_.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
chip.positions<-paste("R0", 1:8, "C01", sep = "")


outOrder<-NULL
IDs<-sample(unique(sampleSheet$Individual))
types<-unique(sampleSheet$Cell.type)
for(each in IDs){
  types<-sample(types)
  for(entry in types){
    selected<-which(sampleSheet$Cell.type == entry & sampleSheet$Individual == each)
    selected<-sample(selected)
    outOrder<-c(outOrder, selected)
  }
}

sampleSheet<-sampleSheet[outOrder,]
sampleSheet$Chip.Location<-chip.positions

write.table(sampleSheet, "SCZHydroxymethFACsSortedPlateDesign14052019.txt")

