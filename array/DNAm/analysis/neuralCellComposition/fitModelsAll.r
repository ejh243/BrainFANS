##---------------------------------------------------------------------#
##
## Title: Derive models for cell composition in brain
##
## Purpose of script: To train a series of models to predict different combinations of neural cell types
##
## Author: Eilis Hannon
##
## Date Created: 19/10/2022
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# identifies CpGs for brain deconvolution models
# uses IDOL method & ANOVA (pickCompProbes)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

set.seed(100)
numProbes<-100 # number of sites included in model
intervals<-seq(0.1,0.9,0.1) # prop of each cell type for test data
probeSelect<- "any" 

args<-commandArgs(trailingOnly = TRUE)
normData<-args[1]
refPath<-args[2]
refPanelPath<-args[3]


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(CETYGO)
library(IDOL)
library(doParallel)
nworkers <- 2 ## if pushed too high causes OOM errors


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#
message("Loading data")
load(normData)
load(file.path(refPath, "AllProbeIlluminaAnno.Rdata"))

colnames(pheno.all)[3]<-"CellType" # need to rename for use with IDOL functions
pheno.all$CellType<-gsub("\\+", "Pos", pheno.all$CellType) ## need to remove the "+" and "-"
pheno.all$CellType<-gsub("\\-", "Neg", pheno.all$CellType) ## need to remove the "+" and "-"
pheno.all$CellType<-factor(pheno.all$CellType)

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


#----------------------------------------------------------------------#
# SELECT PROBES FOR DECONVOLUTION
#----------------------------------------------------------------------#

refPanels<-read.csv(refPanelPath, stringsAsFactors = FALSE)

indexCells <- split(1:nrow(pheno.all), pheno.all$CellType)
exclude<-pheno.all$Basename[unlist(lapply(indexCells, sample, size = 1))]
names(exclude)<-sort(unique(pheno.all$CellType))

brainCoefIDOL<-list()
brainCoefANOVA<-list()
for(i in 1:nrow(refPanels)){
  print(paste0("Model number ", i))
	cellTypes <- unlist(strsplit(refPanels[i,2], ";"))
	cellTypes<-sort(cellTypes)
	pheno.sub<-pheno.all[pheno.all$CellType %in% cellTypes,]
	pheno.sub$CellType<-factor(pheno.sub$CellType)
	norm.sub<-norm.all[,pheno.sub$Basename]
	cellInd<-pheno.sub$CellType
	
	# generate test "bulk profiles"
	matrixSimProp<-expand.grid(rep(list(intervals), length(cellTypes)))
	matrixSimProp<-matrixSimProp[which(signif(rowSums(matrixSimProp),3) == 1),]
	colnames(matrixSimProp)<-cellTypes
	hetBetas <-createBulkProfiles(norm.sub[,exclude[cellTypes]], matrixSimProp) 

	# perform idol dmr test
	idolDMRs<-CandidateDMRFinder.v2(cellTypes, norm.sub, pheno.sub, M = 150,
	   equal.variance = F)		
	## use IDOL method
	if(length(cellTypes) > 2) {
		idolLib<-IDOLoptimize(idolDMRs, hetBetas, matrixSimProp*100, libSize = numProbes*length(cellTypes), numCores = nworkers, maxIt = 300)
		brainCoefIDOL[[i]]<-idolLib$`IDOL Optimized CoefEsts`
	}
	
	## use ANOVA
	compData <- pickCompProbesMatrix(rawbetas=norm.sub, cellInd=cellInd, cellTypes = cellTypes, numProbes = numProbes/2, probeSelect = probeSelect)
	brainCoefANOVA[[i]]<-compData$coefEsts
	save(brainCoefANOVA, brainCoefIDOL, file = file.path(gsub("3_normalised", "4_analysis", dirname(normData)), "CoefBrainModels.rdata"))
}

