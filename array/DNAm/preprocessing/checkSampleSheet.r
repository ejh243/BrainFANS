##---------------------------------------------------------------------#
##
## Title: Check SampleSheet columns
##
## Purpose of script: checks that required sample sheet columns are present 
##                    and preformatted correctly prior to DNAm QC
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# project folder is provided on command line


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

print("checking required sampleSheet colunmns are present and correctly formatted...")

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
dataDir <- "/Users/EW367/testDATA/epicV2sortedBrain/"

configFile <- file.path(dataDir, "/config.r")
source(configFile)

# Column names to test
req_cols <- c('Sample_ID','Individual_ID','Basename', 'test') 					# required column names
opt_cols <- c('Age')																	# optional column names
cond_cols <- c('Sex','Genotype_IID','Cell_Type')	


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(stringdist) # for amatch()

'%ni%' <- Negate('%in%') # define '%ni%' (not in)


#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#


# load sample sheet
sampleSheet<-read.csv(file.path(dataDir,"0_metadata/sampleSheet.csv"), na.strings = c("", "NA"), stringsAsFactors = FALSE)


#----------------------------------------------------------------------#
# CHECK REQUIRED COLUMNS
#----------------------------------------------------------------------#

bad_parameter_exists <- FALSE


check_columns <- function(columns, type_check, warning_message) {
  for (column in columns) {
    
    if(column %ni% colnames(sampleSheet)){
      bad_column_exists <- TRUE
      message("\n'", column, "' column must be in the sampleSheet\n")
      
      # check for close matches
      near <- column[amatch(colnames(sampleSheet), column, maxDist = 2)] # maxDist = n char diff
      
      if(length(near[!is.na(near)])>0){
        nrmatch <- print(data.frame(Expected.colname=near[!is.na(near)], Received.colname=colnames(sampleSheet)[!is.na(near)]))}
      
      } else {
      
  if (!type_check(sampleSheet[[column]])){
    bad_parameter_exists <- TRUE
    message("\n'", column, "' ", warning_message, "\n")
  }
  message("'", column, "' is correctly defined.")
}

  }}


# check required columns exist
check_columns(req_cols, is.character, "must be a string")





if (!exists("ctCheck")) ctCheck <- FALSE
if (ctCheck) {
  message("\n`ctCheck` is set to TRUE. Checking for 'Cell_Type' column")
  check_columns("CellType", is.character, "must be a string")
}

if (!exists("snpCheck")) snpCheck <- FALSE
if (snpCheck) {
  message("\n`snpCheck` is set to TRUE. Cchecking for 'Genotype_IID' column")
  check_columns("Genotype_IID", is.character, "must be a string")
}

if (!exists("sexCheck")) sexCheck <- FALSE
if (sexCheck) {
  message("\n`sexCheck` is set to TRUE. Cchecking for 'Genotype_IID' column")
  check_columns("CellType", is.character, "must be a string")
}

