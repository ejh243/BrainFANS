##---------------------------------------------------------------------#
##
## Title: Profile accuracy of deconvolution of cell composition heterogeneity from simulations
##
## Purpose of script: Trains and tests prediction of specified combination of neural cell types against reconstructed bulk profiles
##
## Author: Eilis Hannon
##
## Date Created: 04/07/2022
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# tests different models for deconvoluting brain cell types
# compares two methods for selecting CpGs to base the deconvolutions on
# compares reference panels consisting of different panels of cell fractions
# compares different numbers of sites for prediction
# where these reference panels are defined in accompanying csv file
# uses houseman CP method to generate cell composition estimates
# test against reconstructed blood profiles where proportions are known


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

set.seed(100)
probeSelect<- "any" 
intervals<-seq(0.1,0.9,0.1) # prop of each cell type for test data
numProbesIter<-seq(20,200, 20) # number of sites included in model
nRepeats<-10

args<-commandArgs(trailingOnly = TRUE)
normData<-args[1]
refPath<-args[2]
refPanelPath<-args[3]
i<-args[4]

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(CETYGO)
library(IDOL)
library(doParallel)
nworkers <- 8 ## if pushed too high causes OOM errors


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#
message("Loading data")
load(normData)
load(file.path(refPath, "AllProbeIlluminaAnno.Rdata"))

colnames(pheno.all)[3]<-"CellType" # need to rename for use with IDOL functions
pheno.all$CellType<-sub("/", "_", pheno.all$CellType) # need to remove "/" from cell type labels

probeAnnot<-probeAnnot[rownames(norm.all),]

## filter to autosome probes
message("Filtering data")
norm.all<-norm.all[!probeAnnot$CHR %in% c("X", "Y"),]

## filter out cross hyb & SNP probes

crosshyb<-read.csv(file.path(refPath, "CrossHybridisingProbesPriceORWeksberg.csv"), row.names = 1)
probes<-read.csv(file.path(refPath, "SNPsinProbesAnno.csv"), row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(norm.all))
remove<-remove[which(is.na(remove) != TRUE)]
norm.all<-norm.all[-remove,]
## remove SNP probes
probes<-probes[row.names(norm.all),]
norm.all<-norm.all[which(probes$Weksburg_CommonSNP_Af_within10bpSBE == "" & probes$Illumina_CommonSNP_Af_within10bpSBE == ""),]

## select reference panel for training
refPanels<-read.csv(refPanelPath, stringsAsFactors = FALSE)
modelNum<-refPanels[i,1]
cellTypes <- unlist(strsplit(refPanels[i,2], ";"))
cellTypes<-sub("/", "_", sort(cellTypes))
pheno.all<-pheno.all[pheno.all$CellType %in% cellTypes,]

norm.all<-norm.all[,pheno.all$Basename]
pheno.all$CellType<-factor(pheno.all$CellType)

indexCells <- split(1:nrow(pheno.all), pheno.all$CellType)
predOut<-NULL

for(j in 1:nRepeats){
	message(paste("Running simulation", j))
	## select one sample to drop from each cell type
	exclude<-unlist(lapply(indexCells, sample, size = 1))
	
	## select one sample from training data to generate bulk training data
	trainIndex<-exclude
	while(sum(trainIndex %in% exclude) > 0){
		trainIndex<-unlist(lapply(indexCells, sample, size = 1))
	}

	pheno.sub<-pheno.all[-exclude,]
	norm.sub<-norm.all[,-exclude]
	cellInd<-pheno.sub$CellType

	# generate training and test "bulk profiles"
	## simulate full spectrum of brain profiles where each cell type is present at at least 10%
	matrixSimProp<-expand.grid(rep(list(intervals), length(cellTypes)))
	matrixSimProp<-matrixSimProp[which(signif(rowSums(matrixSimProp),3) == 1),]
	colnames(matrixSimProp)<-cellTypes

	hetBetas <-createBulkProfiles(norm.all[,trainIndex[cellTypes]], matrixSimProp) ## training
	testBulkBetas <-createBulkProfiles(norm.all[,exclude[cellTypes]], matrixSimProp)		
	# perform idol dmr test with M large enough to only do it once
	idolDMRs<-CandidateDMRFinder.v2(cellTypes, norm.sub, pheno.sub, M = 150,
	   equal.variance = F)		
	message("Selected IDOL DMRs")   
	## test out each algorithm with increasing numbers of sites
	for(numProbes in numProbesIter){
		## select probes as basis of algorithm
		## use ANOVA method
		compData <- pickCompProbesMatrix(rawbetas=norm.sub, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes/2, probeSelect = probeSelect)
		
		## test model in simulated bulk profiles
		predCC<-projectCellTypeWithError(testBulkBetas[rownames(compData$coefEsts),], compData$coefEsts)
		diffCC<-predCC - matrixSimProp
		rmseCC<- sqrt(rowMeans(diffCC^2))
		predOut <- rbind(predOut, cbind("Selection" = "ANOVA", "nProbes" = numProbes*length(cellTypes), matrixSimProp, predCC, diffCC, "RMSE" = rmseCC))
				
		## use IDOL method
		## only works if > 2 cell types
		if(length(cellTypes) > 2) {
			idolLib<-IDOLoptimize(idolDMRs, hetBetas, matrixSimProp*100, libSize = numProbes*length(cellTypes), numCores = nworkers, maxIt = 300)
			predCC<-projectCellTypeWithError(testBulkBetas[rownames(idolLib$`IDOL Optimized CoefEsts`),], idolLib$`IDOL Optimized CoefEsts`)
			diffCC<-predCC - matrixSimProp
			rmseCC<- sqrt(rowMeans(diffCC^2))
			predOut <- rbind(predOut, cbind("Selection" = "IDOL", "nProbes" = numProbes*length(cellTypes), matrixSimProp, predCC, diffCC, "RMSE" = rmseCC))
		}
	}
	save(predOut, file = file.path(gsub("3_normalised", "4_analysis", dirname(normData)), paste0("CompareDeconvolutionParameters_Model", i, ".rdata")))
}

