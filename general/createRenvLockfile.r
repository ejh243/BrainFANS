##---------------------------------------------------------------------#
##
## Title: Create renv files for managing R packages
##
## Purpose of script: To initialise or update the renv files for a 
## set of R scripts in a folder
##
##---------------------------------------------------------------------#

## NOTE this script has to be run interactively and you need to be in the 
## relevant scripts folder

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

install.packages("renv")
library(renv)

#----------------------------------------------------------------------#
# CREATE ENVIRONMENT
#----------------------------------------------------------------------#

renv::init(bare = TRUE)
renv::install("tjgorrie/bigmelon")
renv::install("ds420/CETYGO")
renv::install("EpigeneticsExeter/cdegUtilities")
renv::snapshot()

#----------------------------------------------------------------------#
# UPDATE ENVIRONMENT
#----------------------------------------------------------------------#
renv::snapshot()

