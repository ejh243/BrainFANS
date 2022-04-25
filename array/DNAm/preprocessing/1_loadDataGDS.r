## Written by Eilis and Tyler
## R script to load data from idats files into a gds file
## project folder is provided on command line at execution
## requires a sample sheet called sampleSheet.csv in the 0_metadata folder with Basename or Chip.ID and Chip.Location to identify which samples to load
## assumes idats are in the 1_raw folder
## gds file will be outputted to folder 2_gds
## NB if gds file already exists option to delete and recreate

args<-commandArgs(trailingOnly = TRUE)

library(bigmelon)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(devtools)
devtools::load_all(path = "../functionsR")

dataDir <- args[1]
gdsFile <-paste0(dataDir, "2_gds/raw.gds")
setwd(dataDir) 

## load sample sheet
sampleSheet<-read.csv("0_metadata/sampleSheet.csv", na.strings = c("", "NA"), stringsAsFactors = FALSE)
## if no column Basename, creates from columns Chip.ID and Chip.Location
if(!"Basename" %in% colnames(sampleSheet)){
	sampleSheet$Basename<-paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
}

print(paste(nrow(sampleSheet), "samples identified from sample sheet to be loaded"))


## issue with different numbers of probes on different versions of array
## first classify samples by probe version
nProbes <- sapply(paste0("1_raw/", sampleSheet$Basename, "_Red.idat"), readIDAT, what = "nSNPsRead")
scanDate <- unlist(sapply(paste0("1_raw/", sampleSheet$Basename, "_Red.idat"), getScanDate))

sampleSheet<-cbind(sampleSheet, nProbes, scanDate)

## load data separately
loadGroups<-split(gsub("1_raw/|_Red.idat", "", names(nProbes)), as.factor(nProbes))

for(i in 1:length(loadGroups)){
	gdsFile.sub<-gsub("\\.gds", paste0("_", names(loadGroups)[i], "\\.gds"),gdsFile)
	
	if(file.exists(gdsFile.sub)){
		print(paste("Loading gds file:", gdsFile.sub))
		gfile <- openfn.gds(gdsFile.sub, readonly = FALSE)
		## check which samples already loaded
		sampToLoad <- loadGroups[[i]][!loadGroups[[i]] %in% colnames(gfile)]
		
		if(length(sampToLoad) > 0){
			print(paste("Loading", length(sampToLoad[[i]]), "samples"))
			setwd("1_raw/")
			for(exprID in sampToLoad){
				gfile <- iadd(exprID, gds=gdsFile.sub)
			}
			newSamples<-TRUE
			
		} else {
			print("No new samples to upload")
			newSamples<-FALSE
		}
		updateProbes<-FALSE
	} else {
		print(paste("Creating new gds file", gdsFile.sub))
		print(paste("Loading", length(loadGroups[[i]]), "samples"))
		## load each file and add to gds file
		setwd("1_raw/")
		for(exprID in loadGroups[[i]]){
			gfile <- iadd(exprID, gds=gdsFile.sub)
		}
		updateProbes<-TRUE
		newSamples<-TRUE
	}	

	## update feature data
	if(updateProbes){
		print("Updating Feature data")

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
	}
	
	## update pData need to check all in the sample order, should be as loaded from sample sheet
	## format sampleSheet to match order of gds file
	if(newSamples){
		sampleOrder<-colnames(gfile)
		## check all samples in gfile in samplesheet
		sampleSheet.sub<-sampleSheet[match(sampleOrder, sampleSheet$Basename),]

		sampleSheet.sub<-cbind(sampleSheet.sub$Basename, sampleSheet.sub)
		colnames(sampleSheet.sub)[1]<-"barcode"
		add.gdsn(gfile, 'pData', val = data.frame(lapply(as.data.frame(sampleSheet.sub), as.character), stringsAsFactors = FALSE), replace = TRUE)
	}

	## create back up
	setwd(dataDir)
	## check if any changes
	if(newSamples){
		print("Creating backup of gds file...")
		 f <- createfn.gds(gsub("\\.gds", "_backup.gds", gdsFile.sub))
		 for(node in ls.gdsn(gfile)){
			copyto.gdsn(node = f, source = index.gdsn(gfile, node), name = node)
		}
		closefn.gds(f)
	}
	
	## need to close gds file in order to open in another R session
	closefn.gds(gfile)
}

## merge gds into single files where rownames have been matched
# close any gds that are open
showfile.gds(closeall=TRUE, verbose=TRUE)

mergedf<-createfn.gds(gdsFile)

gfileList<-lapply(paste0(gsub("\\.gds", "" ,gdsFile), "_", names(loadGroups), ".gds"), openfn.gds, readonly = FALSE, allow.fork = TRUE)

## start with betas as function to add rownames and colnames
listBetas <- lapply(lapply(gfileList, betas), read.gdsn)
listProbeIDs <- lapply(lapply(gfileList, index.gdsn, path = "fData/Probe_ID"), read.gdsn)
for(i in 1:length(listBetas)){
	rownames(listBetas[[i]])<-listProbeIDs[[i]]
}

## get list of merged probeids
probeIDs <- unique(unlist(listProbeIDs))

listBetas<-lapply(listBetas, reformatBetas, probeIDs)

allBetas<-do.call(cbind, listBetas)

add.gdsn(mergedf, 'betas', val = allBetas)

## need to merge columns but check rownames in process
for(node in c("pvals","methylated","unmethylated","NBeads")){
	allDat<- lapply(lapply(gfileList, index.gdsn, path = node), read.gdsn)
	for(i in 1:length(allDat)){
		rownames(allDat[[i]])<-listProbeIDs[[i]]
	}
	allDat<-lapply(allDat, reformatBetas, probeIDs)
	allDat <- do.call(cbind, allDat)
	
	add.gdsn(mergedf, node, val = allDat)
}

## add feature data
allFeature <- lapply(lapply(gfileList, index.gdsn, path = "fData"), read.gdsn)

allFeature <- unique(do.call(rbind, allFeature))

## match probe order
allFeature<-allFeature[match(probeIDs, allFeature$Probe_ID),]
add.gdsn(mergedf, "fData", val = allFeature, replace = TRUE)

## add sample data
allPheno <- lapply(lapply(gfileList, index.gdsn, path = "pData"), read.gdsn)

allPheno<-unique(do.call(rbind, allPheno))
add.gdsn(mergedf, "pData", val = allPheno, replace = TRUE)

## don't need to check rownames of these objects
for(node in c("QCmethylated","QCunmethylated")){
	allDat<- lapply(lapply(gfileList, index.gdsn, path = node), read.gdsn)
	allDat <- do.call(cbind, allDat)
	
	add.gdsn(mergedf, node, val = allDat)
}

add.gdsn(mergedf, "QCrownames", val = read.gdsn(index.gdsn(gfileList[[1]], "QCrownames")))

## update history with this merge
add.gdsn(mergedf, "history", val = data.frame("submitted" = Sys.time(), "finished" = Sys.time(), "command" = paste0("Merged", length(listBetas), " gdsfiles with dimensions (", paste(unlist(lapply(lapply(listBetas, dim), paste, collapse=",")), collapse = "),("), ").")))

## set paths to find sample and row names
add.gdsn(mergedf, "paths", val = read.gdsn(index.gdsn(gfileList[[1]], "paths")))

closefn.gds(mergedf)

lapply(gfileList, closefn.gds)

