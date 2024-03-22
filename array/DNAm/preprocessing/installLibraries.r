##INSTALL LIBRARIES SCRIPT
##Installs all R packages required DNAm QC pipeline
#Add in additional lines for Alice's errors... Git and config sourcing?
args<-commandArgs(trailingOnly = TRUE)

source(args[1])

#Bioconductor is required for most installations
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	
	
#Bigmelon required for all scripts
#BiocManager::install("bigmelon")
remotes::install_github("schalkwyk/wateRmelon")
remotes::install_github("tjgorrie/bigmelon")

#GDS script

if(arrayType=='450K'){
  #library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  BiocManager::install(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  BiocManager::install(IlluminaHumanMethylation450kmanifest)
}
if(arrayType=='EPICv1'){
  BiocManager::install(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  BiocManager::install(IlluminaHumanMethylationEPICmanifest)
}
if(arrayType=='EPICv2'){
	install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
	install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
}

install.packages("devtools")

#Calculating QC metrics
install.packages("e1071")


#Creating QC reports
install.packages("pander") 
install.packages("kableExtra")


#Additional packages for Brain Cell Proportion Prediction
BiocManager::install("genefilter")
install.packages("quadprog")


#Packages for CETYGO script
BiocManager::install(c("genefilter", "minfi"))
install.packages("quadprog")

# install devtools to install from GitHub
install.packages("devtools")
library(devtools)

install_github("ds420/CETYGO")

