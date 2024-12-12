
##---------------------------------------------------------------------#
##
## Title: Normalisation
##
## Purpose of script: Perform normalisation of cell-specific DNAm data
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# project folder is provided on command line at execution
# normalisation is performed across whole sample and within cell type
# probe type info from R package was incomplete so loaded separately

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refDir <- args[2]
gdsFile <-file.path(dataDir, "/2_gds/raw.gds")
normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)
qcOutFolder<-file.path(dataDir, "/2_gds/QCmetrics")
normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
configFile <- paste0(dataDir, "/config.r")

source(configFile)

arrayType <- toupper(arrayType)

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(bigmelon, warn.conflicts = FALSE, quietly = TRUE)
library(data.table, warn.conflicts = FALSE, quietly = TRUE)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

probeFilt <- with(new.env(), {load(file.path(qcOutFolder, "QCmetrics.rdata")); get("probeFilt")})


# filter samples
if(ctCheck){
	QCSum<-read.csv(file.path(dataDir, "/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv"), row.names = 1, stringsAsFactors = FALSE)
	passQC<-QCSum$Basename[which(QCSum$passQCS3)]

	QCmetrics<-read.csv(paste0(qcOutFolder,"/QCMetricsPostCellTypeClustering.csv"), stringsAsFactors = FALSE)
	QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]
} else {
	QCSum<-read.csv(file.path(dataDir, "/2_gds/QCmetrics/PassQCStatusAllSamples.csv"), row.names = 1, stringsAsFactors = FALSE)
	passQC<-QCSum$Basename[which(QCSum$passQCS2)]

	QCmetrics<-read.csv(paste0(qcOutFolder,"/QCMetricsPostSampleCheck.csv"), stringsAsFactors = FALSE)
	QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

}
cellTypes<-unique(QCmetrics$Cell_Type)

# create new gfile with only samples that pass QC
if(exists(normgdsFile)){
	file.remove(normgdsFile) ## delete if already exists
}
normfile <- createfn.gds(normgdsFile)
for(node in ls.gdsn(gfile)){
	copyto.gdsn(node = normfile, source = index.gdsn(gfile, node), name = node)
}

# filter out samples that fail QC
# match to basename not colnames
rawbetas<-betas(normfile)[,]
passedProbes <- which(probeFilt[,"pFiltProbesPass"] & probeFilt[,"beadFiltPass"])

subSet(normfile, i=passedProbes, j=match(passQC, colnames(rawbetas)))

# close gds file in order to open in another R session
closefn.gds(gfile)

# need to extract below to run normalisation
meth<-methylated(normfile)[,]
unmeth<-unmethylated(normfile)[,]
rawbetas<-betas(normfile)[,]

manifest <- cdegUtilities::readManifest(
	referenceDirectory = refDir,
	probeMatchingIndex = rownames(rawbetas),
	arrayType = arrayType 
)
if (!exists("manifest"))
	stop("Manifest file could not be loaded correctly")

#----------------------------------------------------------------------#
# NORMALISE ALL SAMPLES TOGETHER
#----------------------------------------------------------------------#

normbeta<-adjustedDasen(
                       onetwo = manifest$designType,
                       chr = manifest$CHR,
                       mns = meth,
                       uns = unmeth)
add.gdsn(normfile, 'normbeta', val = normbeta, replace = TRUE)

#----------------------------------------------------------------------#
# NORMALISE CELL TYPES SEPARATELY
#----------------------------------------------------------------------#

if(length(cellTypes) > 1){

	celltypeNormbeta<-matrix(NA, nrow = nrow(meth), ncol = ncol(meth))
	rownames(celltypeNormbeta)<-rownames(rawbetas)
	colnames(celltypeNormbeta)<-colnames(rawbetas)
	for(each in cellTypes){
		index<-which(QCmetrics$Cell_Type == each)
		if(length(index) > 2){
			celltypeNormbeta[,index] <- as.matrix(adjustedDasen(
                               		onetwo = manifest$designType,
                               		chr = manifest$CHR,
                               		mns = meth[,index],
                               		uns = unmeth[,index]))		
		}
	}
	add.gdsn(normfile, 'celltypenormbeta', val = celltypeNormbeta, replace = TRUE)

	save(celltypeNormbeta, QCmetrics, file = normData)
} else{

	save(normbeta, QCmetrics, file = normData)
}

#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

print(paste0("The final normalised dataset contains ", ncol(rawbetas), " samples and ", nrow(rawbetas), " probes"))

closefn.gds(normfile)

