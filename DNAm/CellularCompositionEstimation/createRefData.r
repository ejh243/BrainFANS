## script to develop neural cell type deconvolution algorithm
## uses houseman method with new ref data for purified brain cell populations
## editted original code from minfi to take a matrix rather than RGset - this means data is not preprocessed together.


source("FunctionsForBrainCellProportionsPrediction.r")
library(minfi)
library(genefilter)
library(rafalib)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
source("rmdConfig.bdr")
setwd(dataDir) 
load(normData)

pheno$Cell.type<-gsub(" \\+ve", "", pheno$Cell.type) ## need to remove the "+" and "-"
pheno$Cell.type<-gsub(" -ve", "", pheno$Cell.type) ## need to remove the "+" and "-"

## parameters
probeSelect = "any" ## options "both" (select equal number of probes associated with hyper and hypo methylation) or "any" (just select probes based on significance regardless of direction)
cellTypes = c("Double","NeuN", "Sox10", "IRF8")
cellTypes<-sort(cellTypes)
numProbes<-100

## select probes as basis of algorithm
compData <- pickCompProbes(rawbetas=celltypenormbeta, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes, probeSelect = probeSelect)
braincelldata <- compData$coefEsts
save(braincelldata, file = "RefDataForCellCompEstimation.rdata")


## re run to select coef probes with 450K data
annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))
newfData<-newfData[rownames(celltypenormbeta),]

compData <- pickCompProbes(rawbetas=celltypenormbeta[which(newfData$Methyl450_Loci == "TRUE"),], cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes, probeSelect = probeSelect)

braincelldata450K<-compData$coefEsts
save(braincelldata450K, file = "RefDataForCellCompEstimation450K.rdata")
