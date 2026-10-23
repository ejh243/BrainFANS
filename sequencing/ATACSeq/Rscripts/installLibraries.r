## ================================================================================##
##                     PIPPA STEP 0.2: Set up R environment                        ##
## ================================================================================##
## EXECUTION: Rscript ./Rscripts/installLibraries.r  <project_config>              ||
## - execute from scripts directory                                                ||
##                                                                                 ||
## DESCRIPTION: This script checks if required R packages are installed and install||
## them if not                                                                     ||
##                                                                                 ||
## INPUTS:                                                                         ||
## - <project_config> : path to config.r file.                                     ||
##                                                                                 ||
## REQUIRES:                                                                       ||
## - R version > 4.3 in conda environment                                          ||
## ================================================================================##

args <- commandArgs(trailingOnly=TRUE)

##Vectors with all required R packages in the pipeline
packagesBase <- c("diptest","segmented")

## Then we install the rest of packages if not installed already
install.packages(setdiff(packagesBase, rownames(installed.packages())))  
install.packages(paste0(configDir,"/config/ptest_1.0-8.tar.gz"))

