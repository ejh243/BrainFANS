##---------------------------------------------------------------------#
##
## Title: Remove Sex Specific Probes
##
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# project folder is provided on command line
# ref directory is provided on command line
# manifest file for epic data is already available in the gds object


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refDir <- args[2]


gdsFile <-paste0(dataDir, "/2_gds/raw.gds")
qcOutFolder<-file.path(dataDir, "/2_gds/QCmetrics")
normgdsFile<-sub("\\.gds", "Norm.gds", gdsFile)
autoData <-paste0(dataDir, "/3_normalised/autoOnlyNormBeta.rdata")
configFile <- paste0(dataDir, "/config.r")
epic2Manifest <- paste0(refDir,"/EPICArray/EPIC-8v2-0_A1.csv")


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

normfile<-openfn.gds(normgdsFile, readonly = FALSE)

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


#----------------------------------------------------------------------#
# LOAD MANIFEST FILE
#----------------------------------------------------------------------#

if(arrayType == "V2"){
manifest<-fread(epic2Manifest, skip=7, fill=TRUE, data.table=F)
manifest<-manifest[match(fData(normfile)$Probe_ID, manifest$IlmnID), c("CHR", "Infinium_Design_Type")]
print("loaded EpicV2 manifest")
}

if(arrayType == "450K"){
load(file.path(refDir, "450K_reference/AllProbeIlluminaAnno.Rdata"))
manifest<-probeAnnot[match(fData(normfile)$Probe_ID, probeAnnot$ILMNID), c("CHR", "INFINIUM_DESIGN_TYPE")]
colnames(manifest) <- c("CHR", "Infinium_Design_Type")
manifest$CHR <- paste0("chr", manifest$CHR)
print("loaded 450K manifest")
rm(probeAnnot)
}


#----------------------------------------------------------------------#
# FILTER TO AUTOSOMAL PROBES ONLY
#----------------------------------------------------------------------#

if(arrayType == "V2" | arrayType == "450K"){
    auto.probes<-which(manifest$CHR != "chrX" & manifest$CHR != "chrY")
} else {
    auto.probes<-which(fData(normfile)$chr != "chrX" & fData(normfile)$chr != "chrY")
}


autoOnlyNormBeta <- betas(normfile)[,][auto.probes,]


#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

add.gdsn(normfile, 'autoOnlyNormBeta', val = autoOnlyNormBeta, replace = TRUE)

save(autoOnlyNormBeta, QCmetrics, file = autoData)

closefn.gds(normfile)
