##---------------------------------------------------------------------#
##
## Title:dmGsea results filtering
##
##
## Input files: results files from running dmGseaAnalysis.R script
##
## Output files: filtered results files
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

dataDir <- "Schizophrenia"
outDir <- file.path(dataDir, "4_analysis/dmGseaV3/") # main results
#outDir <- file.path(dataDir, "4_analysis/dmGseaHypoHyperAnalysis")
#outDir <- file.path(dataDir, "4_analysis/dmGseaNewThresholds")


#----------------------------------------------------------------------#
# LOAD AND FORMAT DATA
#----------------------------------------------------------------------#

library(dplyr)


#----------------------------------------------------------------------#
# function to filter results
#----------------------------------------------------------------------#


# loop through all the results files and count number of pathways with FDR < 0.05

filtRes <- function(runType){
  
  setwd(file.path(outDir, runType))
  
  # read in results files
  fileList <- list.files(".", pattern = runType, full.names = T)
  resFiles<- lapply(fileList, function(x) read.csv(x, stringsAsFactors = F))
  names(resFiles) <- tools::file_path_sans_ext(basename(fileList))
  
  
  # subset to jsut sig results
  sigRes <- lapply(resFiles, function(df) {
    df %>%   dplyr::filter(nGene > 1 & FDR < 0.05) %>%
      arrange(FDR) %>%
      mutate(GO_ID = sub("\\|.*$", "", ID))
  })
  
  
  # save filtered results
  for(name in names(sigRes)){
    outFile <- paste0("filteredResults/", name, "_nGeneMoreThan1.csv")
    write.csv(sigRes[[name]], outFile)
  }
  
}


#----------------------------------------------------------------------#
# Run function
#----------------------------------------------------------------------#

filtRes("gsGene")


