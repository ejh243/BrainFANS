
##INSTALL LIBRARIES SCRIPT
##Installs all R packages required DNAm QC pipeline

##---------------------------------------------------------------------#
##
## Title: Install required packages
##
## Purpose of script: installs all packages required for QC pipeline steps.
##
##---------------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

configFile <- paste0(dataDir, "/config.r")
source(configFile)

#---------------------------------------------------------------------#
# INSTALL BIOCONDUCTOR
#---------------------------------------------------------------------#

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#---------------------------------------------------------------------#
# INSTALL BIGMELON
#---------------------------------------------------------------------#	
	
#Bigmelon required for all scripts
remotes::install_github("schalkwyk/wateRmelon")
remotes::install_github("tjgorrie/bigmelon")

#GDS script

if(arrayType=='450K'){
  #library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  BiocManager::install(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  BiocManager::install(IlluminaHumanMethylation450kmanifest)
}
if(arrayType=='V1'){
  BiocManager::install(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  BiocManager::install(IlluminaHumanMethylationEPICmanifest)
}
if(arrayType=='V2'){
	install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
	install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")

}

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR QC METRICS
#---------------------------------------------------------------------#

install.packages(setdiff(c("e1071","stringdist","data.table"), rownames(installed.packages())),repos = "https://cran.r-project.org")

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR BRAIN CELL PROPORTION PREDICTION
#---------------------------------------------------------------------#
#Creating QC reports
install.packages(setdiff(c("pander","kableExtra"), rownames(installed.packages())),repos = "https://cran.r-project.org")

#Additional packages for Brain Cell Proportion Prediction
BiocManager::install(c("genefilter", "minfi"))
install.packages(setdiff("quadprog", rownames(installed.packages())),repos = "https://cran.r-project.org")

# install devtools to install from GitHub
install.packages(setdiff("devtools", rownames(installed.packages())),repos = "https://cran.r-project.org")
library(devtools)
install_github("ds420/CETYGO")
install_github("EpigeneticsExeter/cdegUtilities")
