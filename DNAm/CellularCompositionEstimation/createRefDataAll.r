
source("FunctionsForBrainCellProportionsPrediction.r")
library(minfi)
library(genefilter)
library(rafalib)

## parameters
probeSelect = "any" ## options "both" (select equal number of probes associated with hyper and hypo methylation) or "any" (just select probes based on significance regardless of direction)
cellTypes = c("GABAneurons","GLIA","GLUneurons", "NeuN","Sox10","IRF8")
cellTypes<-sort(cellTypes)
numProbes<-100

gfile.epi<-openfn.gds("gdsFiles/epigabaNorm.gds")
gfile.bdr<-openfn.gds("gdsFiles/bdrNorm.gds")

QCmetrics.epi<-read.gdsn(index.gdsn(gfile.epi, "QCdata"))
QCmetrics.bdr<-read.gdsn(index.gdsn(gfile.bdr, "QCdata"))
pheno.epi<-as.data.frame(QCmetrics.epi[,c("Basename", "Individual_ID", "Cell.type")], stringsAsFactors = FALSE)
pheno.epi$Cell.type<-gsub(" ", "", pheno.epi$Cell.type)
pheno.bdr<-as.data.frame(QCmetrics.bdr[,c("Basename", "Indidivual.ID", "Cell.type")], stringsAsFactors = FALSE)
pheno.bdr$Cell.type<-gsub(" \\+ve", "", pheno.bdr$Cell.type) ## need to remove the "+" and "-"
pheno.bdr$Cell.type<-gsub(" -ve", "", pheno.bdr$Cell.type) ## need to remove the "+" and "-"

colnames(pheno.bdr)<-colnames(pheno.epi)
pheno.all<-rbind(pheno.epi, pheno.bdr)
pheno.all<-pheno.all[ match(colnames(rawbetas.all), pheno.all$Basename),]

pheno.all$Cell.type<-as.factor(pheno.all$Cell.type)


## use normalised within cell types

norm.epi<-read.gdsn(index.gdsn(gfile.epi, "celltypenormbeta"))
rownames(norm.epi)<-rownames(gfile.epi)
colnames(norm.epi)<-colnames(gfile.epi)
norm.bdr<-read.gdsn(index.gdsn(gfile.bdr, "celltypenormbeta"))
rownames(norm.bdr)<-rownames(gfile.bdr)
colnames(norm.bdr)<-colnames(gfile.bdr)

probes<-intersect(rownames(norm.epi), rownames(norm.bdr))

norm.epi<-norm.epi[probes,]
norm.bdr<-norm.bdr[probes,]

norm.all<-cbind(norm.epi, norm.bdr)


cellInd<-pheno.all$Cell.type

## select probes as basis of algorithm
compData <- pickCompProbes(rawbetas=norm.all, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes, probeSelect = probeSelect)
braincelldata <- compData$coefEsts
save(braincelldata, file = "RefDataForCellCompEstimationAll.rdata")

## exclude NeuN
cellTypes = c("GABAneurons","GLIA","GLUneurons","Sox10","IRF8")
## select probes as basis of algorithm
compData <- pickCompProbes(rawbetas=norm.all, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes, probeSelect = probeSelect)
braincelldata <- compData$coefEsts
save(braincelldata, file = "RefDataForCellCompEstimationAllNoNeuN.rdata")
