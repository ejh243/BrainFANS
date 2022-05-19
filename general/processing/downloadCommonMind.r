 library(synapser) 
 library(synapserutils) 
 args=commandArgs(trailingOnly=TRUE)

 synLogin("ejh243",args[2]) 
 
 setwd("/lustre/projects/Research_Project-MRC190311/")
 
toDownload<-read.table("RNASeq/CommonMind/0_metadata/synapseIDs.txt", stringsAsFactors = FALSE)

## batch into chunks

cutPoints<-round(seq(1, nrow(toDownload), length.out = 11))

print(paste("Running batch", args[1])) 

## CommonMind DLPFC
if(args[1] > 0){
	subIndex<-c(cutPoints[as.numeric(args[1])]:cutPoints[as.numeric(args[1])+1])
	for(each in toDownload[subIndex,1]){
		files <- synapserutils::syncFromSynapse(each, path = "RNASeq/CommonMind/DLPFC", ifcollision="keep.local")
	}
 } else {
  files <- synapserutils::syncFromSynapse("syn18134197", path = "RNASeq/CommonMind/DLPFC", ifcollision="keep.local")
 }