 library(synapser) 
 library(synapserutils) 
 args=commandArgs(trailingOnly=TRUE)

 synLogin("ejh243",args[2]) 
 
 setwd("/lustre/projects/Research_Project-MRC190311/")

## CommonMind DLPFC
if(args[3] == 0){
  files <- synapserutils::syncFromSynapse("syn18134196", path = "RNASeq/CommonMind/DLPFC")
 } else {
  files <- synapserutils::syncFromSynapse("syn18134197", path = "RNASeq/CommonMind/DLPFC")
 }