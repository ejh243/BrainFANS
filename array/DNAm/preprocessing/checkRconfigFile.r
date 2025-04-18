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
configFile <- file.path(dataDir, "config.r")

source(configFile)

qcRmdParams <- c("projectTitle", "processedBy")
qcthres <- c("thresBS", "intenThres", "nvThres", "perMiss")
logicalParams <- c("sexCheck", "snpCheck", "ctCheck")
multimodalParams <- c("xMus", "yMus", "xSigmas", "ySigmas")

ctThres <- c("studentThres", "nSDThres")
ctCellParams <- c("predDistinctCT", "neunCT")



#----------------------------------------------------------------------#
# Check config file parameters
#----------------------------------------------------------------------#
bad_parameter_exists <- FALSE

check_parameters <- function(parameters, type_check, warning_message) {
  for (parameter in parameters) {
    if (!exists(parameter)) {
      bad_parameter_exists <- TRUE
      warning("\n'", parameter, "' must be defined in the config file\n")
    }
    if (!type_check(get(parameter))) {
      bad_parameter_exists <- TRUE
      warning("\n'", parameter, "' ", warning_message, "\n")
    }
    message("'", parameter, "' is correctly defined.")
  }
}

check_parameters(qcthres, is.numeric, "must be numeric")
check_parameters(logicalParams, is.logical, "must be TRUE or FALSE")
check_parameters(multimodalParams, is.numeric, "must be numeric")

if (!exists("ctCheck")) ctCheck <- FALSE
if (ctCheck) {
  message(
    "\n`ctCheck` is set to TRUE.",
    "Additional parameters need to be checked...\n"
  )
  check_parameters(ctCellParams, is.character, "must be a string")
  check_parameters(ctThres, is.numeric, "must be numeric")
}


if (!toupper(tissueType) %in% c("BRAIN", "BLOOD")) {
  bad_parameter_exists <- TRUE
  warning("\nUnrecognised tissueType. Must be either 'blood' or 'brain'\n")
}
if (!toupper(arrayType) %in% c("450K", "V1", "V2")) {
  bad_parameter_exists <- TRUE
  warning("\nUnrecognised arrayType. Must be '450K', 'V1' or 'V2'\n")
}

if (!"Cell_Type" %in% projVar && ctCheck) {
  bad_parameter_exists <- TRUE
  warning("\n'Cell_Type' must be included in projVar if 'ctCheck' is true\n")
}

if (!file.exists(manifestFilePath)) {
  bad_parameter_exists <- TRUE
  warning(
    "\n'manifestFilePath': ", manifestFilePath, " does not exist. ",
    "Please check this file path"
  )
}

if (bad_parameter_exists) {
  stop(
    "\nMalformed config file detected.\n",
    "Please fix variables before running the rest of the pipeline\n"
  )
}

print("\nAll config file parameters are present and correctly formatted\n")
