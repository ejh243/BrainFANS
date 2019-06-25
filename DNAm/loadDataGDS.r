## Written by Eilis and Tyler
## R script to load data from idats files into a gds file
## requires a sample sheet with Basename or Chip.ID and Chip.Location to identify which samples to load
## script to load DNAm data from idat files listed in a sample sheet into gds file using bigmelon
## NB if gds file already exists it is deleted and recreated

iadd <- function (bar, gds, n=T, ...){
    ifile <- basename(bar)
    pieces <- strsplit(ifile, "[_.]")
    slide <- sapply(pieces, '[', 1)
    pos <- sapply(pieces, '[', 2)
    bar <- paste(slide, pos, sep = "_")
    mlu <- methylumIDATepic(barcodes = bar, n=n, ...)
    app2gds(mlu, gds)
}



library(bigmelon)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)

setwd(dataDir)

## load data
## check if gds file exists; if it does it will be deleted and recreated
if(file.exists(gdsFile)){
	file.remove(gdsFile)
}

setwd("iDats")
## load each file from sampleSheet and add to gds file
for(each in sampleSheet$Basename){
	gfile <- iadd(bar = each, gds = paste("../", gdsFile, sep = "")) 
}

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
sampleOrder<-pData(gfile)$barcode
sampleSheet<-sampleSheet[match(sampleOrder, sampleSheet$Basename),]

sampleSheet<-cbind(pData(gfile)$barcode, sampleSheet)
colnames(sampleSheet)[1]<-"barcode"
add.gdsn(gfile, 'pData', val = data.frame(lapply(as.data.frame(sampleSheet), as.character), stringsAsFactors = FALSE), replace = TRUE)


## need to close gds file in order to open in another R session
closefn.gds(gfile)

