## Written by Eilis and Tyler
## R script to load data from idats files into a gds file
## requires a sample sheet with Basename or Chip.ID and Chip.Location to identify which samples to load
## script to load DNAm data from idat files listed in a sample sheet into gds file using bigmelon
## NB if gds file already exists option to delete and recreate
#args<-commandArgs(trailingOnly = TRUE)

library(bigmelon)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)

setwd(dataDir)

## load data
sampToLoad<-sampleSheet$Basename

##check if gds file already exists
if(file.exists(gdsFile)){
	print(paste("Loading gds file:", gdsFile))
	gfile <- openfn.gds(gdsFile, readonly = FALSE)
	## check if all samples are present
	if(sum(sampleSheet$Basename %in% colnames(gfile)) != nrow(sampleSheet)){
	 sampToLoad<-sampToLoad[!sampleSheet$Basename %in% colnames(gfile)]
	}
	## if any missing add to existing gds
	if(length(sampToLoad) > 0){
		print(paste(length(sampToLoad), " samples to load"))
		## load each file from sampleSheet and add to gds file
		setwd("iDats")
		for(each in sampToLoad){
			gfile <- iadd(bar = each, gds = paste("../", gdsFile, sep = ""))
		}
	}
} else { ## otherwise create
	print(paste("Creating gds file:", gdsFile))
	setwd("iDats")
	for(each in sampToLoad){
			gfile <- iadd(bar = each, gds = paste("../", gdsFile, sep = ""))
		}
}	

## as I get error with above if different versions of the array with different numbers of probes here is a work around
#mset <- methylumIDATepic(sampToLoad, force=T, n=T)
#gfile <- es2gds(mset, gdsFile) 


## update feature data

annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))
newfData <- newfData[rownames(gfile), ] # SNP probes will be missing, and be NA’d
rownames(newfData) <- rownames(gfile)
# need to change column name of to ProbeID
colnames(newfData)[which(colnames(newfData) == "Name")]<-"Probe_ID"
newfData$Probe_ID<-rownames(newfData)
add.gdsn(gfile, 'fData', val = data.frame(lapply(as.data.frame(newfData), as.character), stringsAsFactors = FALSE), replace = TRUE)

## update pData need to check all in the sample order, should be as loaded from sample sheet
## format sampleSheet to match order of gds file
sampleOrder<-colnames(gfile)
sampleSheet<-sampleSheet[match(sampleOrder, sampleSheet$Basename),]

sampleSheet<-cbind(sampleSheet$Basename, sampleSheet)
colnames(sampleSheet)[1]<-"barcode"
add.gdsn(gfile, 'pData', val = data.frame(lapply(as.data.frame(sampleSheet), as.character), stringsAsFactors = FALSE), replace = TRUE)


setwd(dataDir)
## create back up

 f <- createfn.gds(gsub("\\.gds", "_backup.gds", gdsFile))
 for(node in ls.gdsn(gfile)){
	copyto.gdsn(node = f, source = index.gdsn(gfile, node), name = node)
}

## need to close gds file in order to open in another R session
closefn.gds(gfile)
closefn.gds(f)




