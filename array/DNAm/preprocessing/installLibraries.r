##INSTALL LIBRARIES SCRIPT
##Installs all R packages required DNAm QC pipeline
#Add in additional lines for Alice's errors... Git and config sourcing?

#Bioconductor is required for most installations
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	
	
#Bigmelon required for all scripts
BiocManager::install("bigmelon")

#GDS script
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") #This one is commented out in original script so unsure if necessary
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
install.packages("devtools")

#Calculating QC metrics
install.packages("e1071")


#Additional packages for Brain Cell Proportion Prediction
BiocManager::install("genefilter")
install.packages("quadprog")

