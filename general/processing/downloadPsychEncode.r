library(synapser) 
library(synapserutils) 
args=commandArgs(trailingOnly=TRUE)


synLogin('ejh243', args[1]) 
 
setwd("/lustre/projects/Research_Project-MRC190311/")
 
## BrainGVEX
#files <- synapserutils::syncFromSynapse("syn8113123", path = "ATACSeq/RawData/GVEX") 
#files <- synapserutils::syncFromSynapse("syn7062404", path = "RNASeq/RawData/GVEX", ifcollision="keep.local")
#files <- synapserutils::syncFromSynapse("syn5588763", path = "Phenotype/GVEX")
#files <- synapserutils::syncFromSynapse("syn5592792", path = "Phenotype/GVEX")
#files <- synapserutils::syncFromSynapse("syn7105868", path = "Phenotype/GVEX")
 
## EpiDiff/EpiMap (Akbarian)
#files <- synapserutils::syncFromSynapse("syn5691267", path = "ChipSeq/AlignedData/EpiMap") 
#files <- synapserutils::syncFromSynapse("syn5691266", path = "Phenotype/EpiMap")
  
## EPIGABA
#files <- synapserutils::syncFromSynapse("syn12033254", path = "ChIPSeq/epiGaba/1_raw") 
#files <- synapserutils::syncFromSynapse("syn7072866", path = "DNAm/iDats/epiGaba") 
#files <- synapserutils::syncFromSynapse("syn4588489", path = "Phenotype/epiGaba")  
#files <- synapserutils::syncFromSynapse("syn12033256", path = "WGBS/epiGaba/1_raw") 
#files <- synapserutils::syncFromSynapse("syn17096984", path = "oxBS/epiGaba/1_raw") 
#files <- synapserutils::syncFromSynapse("syn4864491", path = "ERRBS/epiGaba/1_raw")

files <- synGet(args[2], downloadLocation = args[3])