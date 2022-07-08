##---------------------------------------------------------------------#
##
## Title: Normalisation
##
## Purpose of script: Perform normalisation of cell-specific DNAm data
##
## Author: Eilis Hannon
##
## Date Created: 2020
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

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
library(bigmelon)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

# filter samples
QCSum<-read.csv(file.path(dataDir, "/2_gds/QCmetrics/passQCStatusStage3AllSamples.csv"), row.names = 1, stringsAsFactors = FALSE)
passQC<-QCSum$Basename[which(QCSum$passQCS3)]

QCmetrics<-read.csv(paste0(qcOutFolder,"/QCMetricsPostCellTypeClustering.csv"), stringsAsFactors = FALSE)
QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

cellTypes<-unique(QCmetrics$Cell.type)

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
subSet(normfile, i=1:length(rownames(normfile)), j=match(passQC, colnames(rawbetas))) 

# close gds file in order to open in another R session
closefn.gds(gfile)

# need to extract below to run normalisation
meth<-methylated(normfile)[,]
unmeth<-unmethylated(normfile)[,]
rawbetas<-betas(normfile)[,]

# need to know which probe type
probeAnno<-read.table(file.path(refDir, "EPICArray/EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
probeAnno<-probeAnno[match(rownames(rawbetas), probeAnno$probeID),]

#----------------------------------------------------------------------#
# NORMALISE ALL SAMPLES TOGETHER
#----------------------------------------------------------------------#

#dasen(normfile, node="normbeta", onetwo=probeAnno$designType)

normbeta<-dasen(meth, unmeth, probeAnno$designType)
add.gdsn(normfile, 'normbeta', val = normbeta, replace = TRUE)

#----------------------------------------------------------------------#
# NORMALISE CELL TYPES SEPARATELY
#----------------------------------------------------------------------#

celltypeNormbeta<-matrix(NA, nrow = nrow(meth), ncol = ncol(meth))
rownames(celltypeNormbeta)<-rownames(rawbetas)
colnames(celltypeNormbeta)<-colnames(rawbetas)
for(each in cellTypes){
	index<-which(QCmetrics$Cell.type == each)
	if(length(index) > 2){
		celltypeNormbeta[,index]<-dasen(meth[,index], unmeth[,index], probeAnno$designType)
	}
}

#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

add.gdsn(normfile, 'celltypenormbeta', val = celltypeNormbeta, replace = TRUE)

closefn.gds(normfile)

save(celltypeNormbeta, QCmetrics, file = normData)