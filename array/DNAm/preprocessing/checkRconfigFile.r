## ---------------------------------------------------------------------#
##
## Title: Check config.r
##
## Purpose of script: Check config file parameters are present and in correct format
##                    and if not exit the script with error message
##
## ---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# Project Folder is provided on the command line


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

print("checking config.r file parameters are present and correctly formatted...")

args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
configFile <- paste0(dataDir, "/config.r")

source(configFile)


# group config variables by type

qcRmdParams <- c("projectTitle", "processedBy")

qcthres <- c("thresBS", "intenThres", "nvThres", "perMiss")

logicalParams <- c("sexCheck", "snpCheck", "ctCheck")

ungrouped <- c("tissueType", "arrayType", "techVar", "bioVar")

# only used if ctCheck = TRUE
ctThres <- c("studentThres", "nSDThres")

ctCellParams <- c("predDistinctCT", "neunCT")



#----------------------------------------------------------------------#
# Check config file parameters
#----------------------------------------------------------------------#

check_parameters <- function(parameters, type_check, stop_message) {
  for (parameter in parameters) {
    if (!exists(parameter)) {
      stop("'", parameter, "' must be defined in the config file")
    }
    if (!type_check(get(parameter))) {
      stop("'", parameter, "' ", stop_message)
    }
    message("'", parameter, "' is correctly defined.")
  }
}

check_parameters(qcthres, is.numeric, "must be numeric")
check_parameters(logicalParams, is.logical, "must be TRUE or FALSE")

if (ctCheck) {
  message("ctCheck is set to TRUE, additional parameters need to be checked..")
  check_parameters(ctCellParams, is.character, "must be a string")
  check_parameters(ctThres, is.numeric, "must be numeric")
}

# check 'ungrouped' params are correctly formatted

if (!toupper(tissueType) %in% c("BRAIN", "BLOOD")) {
  stop("Unrecognised tissueType. Must be either 'blood' or 'brain'")
}


if (!toupper(arrayType) %in% c("HM450K", "V1", "V2")) {
  stop("Unrecognised arrayType. Must be 'HM450K', 'V1' or 'V2'")
}


for (i in c("Sentrix_ID", "Sentrix_Position")) {
  if (!i %in% techVar) {
    stop("'", i, "' must be included in techVar")
  }
}

for (i in c("Individual_ID", "Cell_Type", "Sex")) {
  if (!i %in% bioVar) {
    stop("'", i, "' must be included in bioVar")
  }
}

for (i in c("Sentrix_ID", "Sentrix_Position")) {
  if (!i %in% techVar) {
    stop("'", i, "' must be included in techVar")
  }
}



print("All config file parameters are present and corerctly formatted")
