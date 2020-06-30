## script to develop neural cell type deconvolution algorithm
## uses houseman method with new ref data for purified brain cell populations
## editted original code from minfi to take a matrix rather than RGset - this means data is not preprocessed together.


source("FunctionsForBrainCellProportionsPrediction.r")
library(minfi)
library(genefilter)
library(rafalib)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
source("rmdConfig.epigaba")
setwd(dataDir) 

gfile.epi<-openfn.gds("gdsFiles/epigabaNorm.gds")
norm.epi<-read.gdsn(index.gdsn(gfile.epi, "celltypenormbeta"))
rownames(norm.epi)<-rownames(gfile.epi)
colnames(norm.epi)<-colnames(gfile.epi)

pheno.epi<-read.gdsn(index.gdsn(gfile.epi, "QCdata"))

pheno.epi$Cell.type<-gsub(" ", "", pheno.epi$Cell.type)

## parameters
probeSelect = "any" ## options "both" (select equal number of probes associated with hyper and hypo methylation) or "any" (just select probes based on significance regardless of direction)
cellTypes = c("GABAneurons","GLIA","GLUneurons")
cellTypes<-sort(cellTypes)
numProbes<-100
cellInd<-pheno.epi$Cell.type

## select probes as basis of algorithm
compData <- pickCompProbes(rawbetas=norm.epi, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes, probeSelect = probeSelect)
braincelldata <- compData$coefEsts
save(braincelldata, file = "RefDataEpiGABAForCellCompEstimation.rdata")

