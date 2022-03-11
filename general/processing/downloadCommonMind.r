 library(synapser) 
 library(synapserutils) 
 args=commandArgs(trailingOnly=TRUE)

 synLogin("ejh243",args[2]) 
 
 setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/")

## CommonMind DLPFC
  files <- synapserutils::syncFromSynapse("syn18134196", path = "RNASeq/commonMind/1_raw/DLPFC")
  files <- synapserutils::syncFromSynapse("syn18358503", path = "ATACSeq/commonMind/1_raw/DLPFC", ifcollision="keep.local") 

 files <- synapserutils::syncFromSynapse("syn3275213", path = "ATACSeq/commonMind/0_metadata") 
syn18134196