##---------------------------------------------------------------------#
##
## Title: Check config.r 
##
## Purpose of script: Check config file parameters are present and in correct format
##                    and if not exit the script with error message
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# Project Folder is provided on the command line


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

print("checking config.r file parameters are present and correctly formatted...")

'%ni%' <- Negate('%in%') # define '%ni%' (not in)

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
configFile <- paste0(dataDir, "/config.r")

source(configFile)


# group config variables by type

qcRmdParams <- c("projectTitle",
                 "processedBy")           # strings



qcthres <- c("thresBS",
              "intenThres",
              "nvThres",
              "perMiss") # numeric


logicalParams <- c("sexCheck",
                   "snpCheck",
                   "ctCheck")             # logical



# ungrouped params
ungrouped <- c("tissueType", "arrayType", "techVar", "bioVar")



# only used if ctCheck = TRUE
ctThres <- c("studentThres",         
             "nSDThres")                  # numeric

ctCellParams <- c("predDistinctCT",               
                  "neunCT")               # strings 



#----------------------------------------------------------------------#
# Check config file parameters
#----------------------------------------------------------------------#


# check necessary parameters exist

for(i in c(qcRmdParams, qcthres, logicalParams, ungrouped)){
  if(!exists(i)){
    stop(paste0("'", i, "' must be defined in the config file" ))
  }
}


# if ctCheck is TRUE, check extra params and thresholds needed to run clusterCellTypes.r script

if(ctCheck == TRUE){
  for(i in c(ctCellParams, ctThres)){
    if(!exists(i)){
      stop(paste0("To run ctCheck'", i, "' must be defined in the config file" ))
    }
  }
}


# check calcQCmetrics.r script thresholds

for(i in qcthres){
    if(!is.numeric(get(i))){
      stop(paste0("'", i, "' must be numeric"))
  } 
} 



# check logical parameters

for(i in logicalParams){
    if(!is.logical(get(i))){
      stop(paste0("'", i, "' must be True/False"))
  } 
} 

                    

# if ctCheck is TRUE, check extra params and thresholds needed to run clusterCellTypes.r script


if(ctCheck == TRUE){
  for(i in ctCellParams){
    if(exists(i)){
      if(!is.character(get(i))){
        stop(paste0("'", i, "' must be string"))
      }
    } else {
      stop(paste0("To run ctCheck '", i, "' must be defined in the config file" ))
    }
  }
  for(i in ctThres){
    if(exists(i)){
      if(!is.numeric(get(i))){
        stop(paste0("'", i, "' must be numeric"))
      } 
    } else {
      stop(paste0("To run ctCheck '", i, "' must be defined in the config file" ))
    }
  } 
}


# check 'ungrouped' params are correctly formatted

if(!toupper(tissueType) %in% c("BRAIN", "BLOOD")){
  stop("Unrecognised tissueType. Must be either 'blood' or 'brain'")
}


if(!toupper(arrayType) %in% c("HM450K", "V1", "V2")){
  stop("Unrecognised arrayType. Must be 'HM450K', 'V1' or 'V2'")
}


for(i in c("Sentrix_ID", "Sentrix_Position")){
  if(!i %in% techVar){
    stop(paste0("'", i, "' must be included in techVar"))
  }
}


for(i in c("Individual_ID","Cell_Type","Sex")){
  if(!i %in% bioVar){
    stop(paste0("'", i, "' must be included in bioVar"))
  }
}



print("All config file parameters are present and corerctly formatted")

