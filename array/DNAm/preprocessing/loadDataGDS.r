##---------------------------------------------------------------------#
##
## Title: Load data into GDS object
##
## Purpose of script: from a given sample sheet load data from idats into GDS file
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# project folder is provided on command line at execution
# assumes idats are in 1_raw

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

gdsFile <- file.path(dataDir, "2_gds/raw.gds")

configFile <- paste0(dataDir, "/config.r")
source(configFile)

arrayType <- toupper(arrayType)

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(bigmelon, warn.conflicts = FALSE, quietly = TRUE)
library(cdegUtilities, warn.conflicts = FALSE, quietly = TRUE)


if(arrayType=='450K'){
  library("IlluminaHumanMethylation450kanno.ilmn12.hg19", warn.conflicts = FALSE, quietly = TRUE)
  library("IlluminaHumanMethylation450kmanifest", warn.conflicts = FALSE, quietly = TRUE)
}
if(arrayType=='V1'){
  library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", warn.conflicts = FALSE, quietly = TRUE)
  library("IlluminaHumanMethylationEPICmanifest", warn.conflicts = FALSE, quietly = TRUE)
}
if(arrayType=='V2'){
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38, warn.conflicts = FALSE, quietly = TRUE)
library(IlluminaHumanMethylationEPICv2manifest, warn.conflicts = FALSE, quietly = TRUE)

}

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
setwd(dataDir)
## load sample sheet
sampleSheet <- read.csv("0_metadata/sampleSheet.csv", na.strings = c("", "NA"), stringsAsFactors = FALSE)
## if no column Basename, creates from columns Chip.ID and Chip.Location
if(!"Basename" %in% colnames(sampleSheet)){
  sampleSheet$Basename <- paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
}

print(paste(nrow(sampleSheet), "samples identified from sample sheet to be loaded"))

## check unique experimental IDs

if(nrow(sampleSheet) == length(unique(sampleSheet$Basename))){
  message("All experimental IDs (Basename) are unique")
} else {
  stop("duplicate experimental IDs (Basename) present. Please check sample sheet")
}


## issue with different numbers of probes on different versions of array
## first classify samples by probe version
nProbes <- sapply(paste0("1_raw/", sampleSheet$Basename, "_Red.idat"), readIDAT, what = "nSNPsRead")
if(length(nProbes)==0){
  stop("Error calculating number of probes from IDATs.")
}
sampleSheet <- cbind(sampleSheet, nProbes)

tryCatch(
  expr= {
    scanDate <- unlist(vapply(
      paste0("1_raw/", sampleSheet[["Basename"]], "_Red.idat"),
      cdegUtilities::getScanDate,
      character(1)
    ))
    sampleSheet <- cbind(sampleSheet, scanDate)
  },
  error = function(e) {
    print(e)
    message("No scan date could be found in at least one IDAT file.")
  }
)

## load data separately
loadGroups <- split(gsub("1_raw/|_Red.idat", "", names(nProbes)), as.factor(nProbes))

for(i in 1:length(loadGroups)){
  setwd(dataDir)
  gdsFile.sub <- gsub("\\.gds", paste0("_", names(loadGroups)[i], "\\.gds"), gdsFile)

  ## check if partially-created or fully-created GDS file exists already
  if(file.exists(gdsFile.sub) || file.exists(gdsFile)){

    ## stop if both versions present
    if(file.exists(gdsFile.sub) && file.exists(gdsFile)){
      stop(paste0("Identified 2 versions of GDS file in", dataDir, "/2_gds/. Please (re)move one before proceeding."))
    }

    present <- TRUE

    ## determine which version is present
    if(file.exists(gdsFile)){
      flnm <- gdsFile
    }
    if(file.exists(gdsFile.sub)){
      flnm <- gdsFile.sub
    }

  } else {
      present <- FALSE
      flnm <- gdsFile.sub
  }

  if(present){
    print(paste("Loading gds file:", flnm))
    gfile <- openfn.gds(flnm, readonly = FALSE)
    ## check which samples already loaded
    sampToLoad <- loadGroups[[i]][!loadGroups[[i]] %in% colnames(gfile)]

    if(length(sampToLoad) > 0){
      print(paste("Loading", length(sampToLoad), "samples"))
      setwd("1_raw/")
      for(exprID in sampToLoad){
        gfile <- iadd(exprID, gds=flnm)
      }
      newSamples<-TRUE

    } else {
      print("No new samples to upload")
      newSamples <- FALSE
    }
    updateProbes <- FALSE
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
    if(arrayType== "450K"){
      annoObj <- minfi::getAnnotationObject("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    }
    if(arrayType == "V1"){
      annoObj <- minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    }
    if(arrayType == "V2"){
      annoObj <- minfi::getAnnotationObject("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
    }
    all <- minfi:::.availableAnnotation(annoObj)$defaults
    newfData <- do.call(cbind, lapply(all, function(wh) {
      minfi:::.annoGet(wh, envir = annoObj@data)
    }))
    newfData <- newfData[rownames(gfile), ] # SNP probes will be missing, and be NA’d
    rownames(newfData) <- rownames(gfile)
    # need to change column name of to ProbeID
    colnames(newfData)[which(colnames(newfData) == "Name")] <- "Probe_ID"
    newfData$Probe_ID <- rownames(newfData)
    add.gdsn(gfile, 'fData', val = data.frame(lapply(as.data.frame(newfData), as.character), stringsAsFactors = FALSE), replace = TRUE)
  }

  ## update pData need to check all in the sample order, should be as loaded from sample sheet
  ## format sampleSheet to match order of gds file
  if(newSamples){
    sampleOrder <- colnames(gfile)
    ## check all samples in gfile in samplesheet
    sampleSheet.sub <- sampleSheet[match(sampleOrder, sampleSheet$Basename), ]

    sampleSheet.sub <- cbind(sampleSheet.sub$Basename, sampleSheet.sub)
    colnames(sampleSheet.sub)[1]<-"barcode"
    add.gdsn(gfile, 'pData', val = data.frame(lapply(as.data.frame(sampleSheet.sub), as.character), stringsAsFactors = FALSE), replace = TRUE)
  }
  ## need to close gds file in order to open in another R session
  closefn.gds(gfile)
}
# close any gds that are open
showfile.gds(closeall=TRUE, verbose=TRUE)

#----------------------------------------------------------------------#
# MERGE GDS
#----------------------------------------------------------------------#
if(length(loadGroups) > 1){
  mergedf <- createfn.gds(gdsFile)

  gfileList <- lapply(paste0(gsub("\\.gds", "" ,gdsFile), "_", names(loadGroups), ".gds"), openfn.gds, readonly = FALSE, allow.fork = TRUE)

  ## start with betas as function to add rownames and colnames
  listBetas <- lapply(lapply(gfileList, betas), read.gdsn)
  listProbeIDs <- lapply(lapply(gfileList, index.gdsn, path = "fData/Probe_ID"), read.gdsn)
  listSampleIDs <- lapply(gfileList, colnames)
  for(i in 1:length(listBetas)){
    rownames(listBetas[[i]]) <- listProbeIDs[[i]]
	colnames(listBetas[[i]]) <- listSampleIDs[[i]]
  }

  ## get list of merged probeids
  probeIDs <- unique(unlist(listProbeIDs))

  listBetas <- lapply(listBetas, reformatBetas, probeIDs)

  allBetas <- do.call(cbind, listBetas)

  add.gdsn(mergedf, 'betas', val = allBetas)

  ## need to merge columns but check rownames in process
  for(node in c("pvals","methylated","unmethylated","NBeads")){
    allDat <- lapply(lapply(gfileList, index.gdsn, path = node), read.gdsn)
    for(i in 1:length(allDat)){
      rownames(allDat[[i]]) <- listProbeIDs[[i]]
    }
    allDat <- lapply(allDat, reformatBetas, probeIDs)
    allDat <- do.call(cbind, allDat)

    add.gdsn(mergedf, node, val = allDat)
  }

  ## add feature data
  allFeature <- lapply(lapply(gfileList, index.gdsn, path = "fData"), read.gdsn)

  allFeature <- unique(do.call(rbind, allFeature))

  ## match probe order
  allFeature <- allFeature[match(probeIDs, allFeature$Probe_ID),]
  add.gdsn(mergedf, "fData", val = allFeature, replace = TRUE)

  ## add sample data
  allPheno <- lapply(lapply(gfileList, index.gdsn, path = "pData"), read.gdsn)

  allPheno <- unique(do.call(rbind, allPheno))
  add.gdsn(mergedf, "pData", val = allPheno, replace = TRUE)

  ## don't need to check rownames of these objects
  for(node in c("QCmethylated","QCunmethylated")){
    allDat <- lapply(lapply(gfileList, index.gdsn, path = node), read.gdsn)
    allDat <- do.call(cbind, allDat)

    add.gdsn(mergedf, node, val = allDat)
  }

  add.gdsn(mergedf, "QCrownames", val = read.gdsn(index.gdsn(gfileList[[1]], "QCrownames")))

  ## update history with this merge
  add.gdsn(mergedf, "history", val = data.frame("submitted" = Sys.time(), "finished" = Sys.time(), "command" = paste0("Merged", length(listBetas), " gdsfiles with dimensions (", paste(unlist(lapply(lapply(listBetas, dim), paste, collapse=",")), collapse = "),("), ").")), replace = TRUE)

  ## set paths to find sample and row names
  add.gdsn(mergedf, "paths", val = read.gdsn(index.gdsn(gfileList[[1]], "paths")), replace = TRUE)

  closefn.gds(mergedf)

  lapply(gfileList, closefn.gds)
} else{
  ## if GDS did not originally exist or was partially-created, now rename to raw.gds
  if(flnm == gdsFile.sub){
    print(paste("Renaming file to", gdsFile))
    file.rename(gdsFile.sub, gdsFile)
  }
}
