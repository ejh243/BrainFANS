##---------------------------------------------------------------------#
##
## Title: Install required packages
##
## Purpose of script: installs all packages required for QC pipeline steps.
##
## Author: Rhiannon Haigh
##
## Date Created: 08/2023
##
##---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# INSTALL BIOCONDUCTOR
#---------------------------------------------------------------------#

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#---------------------------------------------------------------------#
# INSTALL BIGMELON
#---------------------------------------------------------------------#	
	
#Bigmelon required for all scripts
BiocManager::install("bigmelon")

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR LOADING GDS FILES
#---------------------------------------------------------------------#

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") #This one is commented out in original script so unsure if necessary
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
install.packages("devtools")

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR QC METRICS
#---------------------------------------------------------------------#

install.packages("e1071")

#---------------------------------------------------------------------#
# INSTALL PACKAGES FOR BRAIN CELL PROPORTION PREDICTION
#---------------------------------------------------------------------#

BiocManager::install("genefilter")
install.packages("quadprog")

