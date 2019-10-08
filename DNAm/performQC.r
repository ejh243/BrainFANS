## this script brings together:
## the script to load the DNAm data from idats into a gds file
## the script that generates QC metrics
## it is written to minimise the need to reload/re perform data QC 
## open and close the gds file within each script

args<-commandArgs(trailingOnly = TRUE)

library(bigmelon)
if (length(args)==0) {
  stop("Config file missing on command line", call.=FALSE)
} else {
  print(paste("Using the following config file:", args[1]))
}

source(args[1])



setwd(dataDir) 

## load sample sheet
sampleSheet<-read.csv(sampleFile, na.strings = "")
## if no column Basename, creates from columns Chip.ID and Chip.Location
if(!"Basename" %in% colnames(sampleSheet)){
	sampleSheet$Basename<-paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
}

print(paste(nrow(sampleSheet), "samples identified from sample sheet to be loaded"))

## to avoid loading the data multiple times check if gds files exist

if(file.exists(gdsFile)){
	gfile <- openfn.gds(gdsFile)
	print(paste("Loading gds file:", gdsFile))


	## does it contain all the samples we are interested in?
	if(sum(sampleSheet$Basename %in% pData(gfile)$barcode) < nrow(sampleSheet)){
		#closefn.gds(gfile)
		source(paste(scriptFolder, "/loadDataGDS.r", sep = "")) ## reload data if some samples are missing
	}
	## close gds file
	closefn.gds(gfile)	
} else { ## otherwise create
	source(paste(scriptFolder, "/loadDataGDS.r", sep = ""))
	print(paste("Creating gds file:", gdsFile))
}


## check if QC metrics have been generated
if(file.exists(qcData)){
	load(qcData)
	
	## check if QC metrics present for all samples
	if(nrow(QCmetrics) != nrow(sampleSheet)){
		file.remove(qcData)
		source(paste(scriptFolder, "/calcQCMetrics.r", sep = ""))
	}
} else {
	source(paste(scriptFolder, "/calcQCMetrics.r", sep = ""))
}



