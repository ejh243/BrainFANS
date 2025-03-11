## Study parameters
projectTitle <- "MRC Schizophrenia FANS samples"
processedBy <- "Complex Disease Epigenetic Group, University of Exeter Medical School"
tissueType <- "brain" # needs to be brain or blood to generate relevant cell composition estimates
arrayType <- "V2" # must be one of "450K" / "V1" / "V2"

## Specify the full path to the manfiest file for your data
manifestFilePath <- ""

## project variables
projVar <- c("Cell_Type", "Sex") # Used to colour plots


## QC thresholds
thresBS <- 80
intenThres <- 500
nvThres <- 0.1
perMiss <- 2

## Multimodality parameters for sex prediction
## To be used with creating mixtures of normal distributions
## Don't change these unless sex prediction incorrectly fails for your dataset
xMus<-c(0.99,1.01) # means
xSigmas<-c(0.05,0.05) # standard deviations
yMus<-c(0.3,1.02)
ySigmas<-c(0.2,0.2)

## Optional QC steps
sexCheck <- TRUE
snpCheck <- TRUE
ctCheck <- TRUE


## ctCheck variables
predDistinctCT <- c("NeuN+", "Sox10+", "IRF8+", "Triple-", "Total")
neunCT <- c("NeuN+", "SATB2+")
studentThres <- 1.5
nSDThres <- 3
