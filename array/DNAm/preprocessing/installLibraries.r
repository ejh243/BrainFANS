
##INSTALL LIBRARIES SCRIPT
##Installs all R packages required DNAm QC pipeline

##---------------------------------------------------------------------#
##
## Title: Install required packages
##
## Purpose of script: installs all packages required for QC pipeline steps.
##
##---------------------------------------------------------------------#

'%ni%' <- Negate('%in%')

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
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  BiocManager::install("IlluminaHumanMethylation450kmanifest")
}
if(arrayType=='V1'){
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
}
if(arrayType=='V2'){
	install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
	install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")

}

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR QC METRICS
#---------------------------------------------------------------------#
pkgs_qc <- c("e1071","stringdist","data.table")
install.packages(setdiff(pkgs_qc, rownames(installed.packages())),repos = "https://cran.r-project.org")

#pipeline relies on matrixStats version 1.1.0
if(packageVersion("matrixStats") != "1.1.0") {
  message("matrixStats version is > 1.1.0, downgrading to version 1.1.0")
  remotes::install_version("matrixStats", version = "1.1.0", quiet = TRUE)
}
packageVersion("matrixStats")

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR BRAIN CELL PROPORTION PREDICTION
#---------------------------------------------------------------------#
#Creating QC reports
pkgs_rep <- c("pander","kableExtra")
install.packages(setdiff(pkgs_rep, rownames(installed.packages())),repos = "https://cran.r-project.org")

#Additional packages for Brain Cell Proportion Prediction
BiocManager::install(c("genefilter", "minfi"))
pkgs_pred <- c("devtools","quadprog")
install.packages(setdiff(pkgs_pred, rownames(installed.packages())),repos = "https://cran.r-project.org")

# load devtools to install from GitHub
library(devtools)
install_github("ds420/CETYGO")
install_github("EpigeneticsExeter/cdegUtilities")


## Check all packages installed successfully ##
all_pkgs <- c("wateRmelon", "bigmelon", pkgs_qc, pkgs_rep, "genefilter", "minfi", pkgs_pred, "CETYGO", "cdegUtilities")
if(arrayType=='450K'){
  all_pkgs <- c(all_pkgs, "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylation450kmanifest")
}
if(arrayType=='V1'){
  all_pkgs <- c(all_pkgs, "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", "IlluminaHumanMethylationEPICmanifest")
}
if(arrayType=='V2'){
  all_pkgs <- c(all_pkgs, "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "IlluminaHumanMethylationEPICv2manifest")
}
if(all(all_pkgs %in% rownames(installed.packages()))){
  print("All packages successfully installed")
}else{
  absent <- all_pkgs[all_pkgs %ni% rownames(installed.packages())]
  print(paste("Failed installation of", length(absent), "packages:", absent))
}




