##---------------------------------------------------------------------#
##
## Title: Format data for MatrixEQTL 
##
## Purpose of script: match DNAm, genotype and covariate data ready
## for matrix QTL analysis
##
## Author: Eilis Hannon
##
## Date Created: 2022-08-11
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refDir<-args[2]
genoDir<-args[3]
pop<-args[4]

normData<-file.path(dataDir, "3_normalised", "normalised.rdata")
resPath<-file.path(dataDir, "4_analysis", "QTLs", "Input", pop)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
cellTypes<-unique(QCmetrics$Cell.type)

load(file.path(refDir, "AFSnpProbesCrossHybProbesToExcludeEPIC.rdata"))

celltypeNormbeta<-celltypeNormbeta[!rownames(celltypeNormbeta) %in% remove,]

probeAnnot<-read.table(file.path(refDir, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE)
probeAnnot<-probeAnnot[which(probeAnnot$chrm != "*"),]
celltypeNormbeta<-celltypeNormbeta[rownames(celltypeNormbeta) %in% probeAnnot$probeID,]
probeAnnot<-probeAnnot[match(rownames(celltypeNormbeta), probeAnnot$probeID),]
probeAnnot$chrm<-gsub("chr", "", as.character(probeAnnot$chrm))
probeAnnot<-probeAnnot[order(probeAnnot$start),]

genoFiles<-list.files(file.path(genoDir, "5_qtls", pop, "hg38"), pattern = "_rsq0.3_[0-9]*.traw", full.names = TRUE)
message("Found the following genotype files")
print(genoFiles)
# read enough of a geno file to match IDs
genoRaw<-read.table(genoFiles[1], header = T, stringsAsFactors = FALSE, nrows = 6)
colIDs<-sub("X", "", unlist(lapply(lapply(strsplit(colnames(genoRaw), "_"), head, n = 2), paste, collapse = "_")))
idMap<-read.table(file.path(genoDir, "0_metadata", "UpdateIDs.txt"), stringsAsFactors = FALSE)


# filter to samples with geno data only
QCmetrics<-QCmetrics[!is.na(match(QCmetrics$Genotype.IID, idMap$V4[match(colIDs, idMap$V3)])),]
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

iIndex<-match(QCmetrics$Genotype.IID, idMap$V4[match(colIDs, idMap$V3)])

message("Processing genotype and methylation data")

for(each in genoFiles){
	#extract chr
	chr <- gsub("\\.traw$", "", tail(unlist(strsplit(basename(each), "_")), n = 1))
	
	message(paste("Running chr", chr))
	# match geno data
	genoRaw<-read.table(each, header = T, stringsAsFactors = FALSE)
	rownames(genoRaw)<-genoRaw$SNP
	genoMap<-genoRaw[,c("SNP", "CHR", "POS")]
	colnames(genoMap)<-c("snp", "chr", "pos")
	genoRaw<-genoRaw[,iIndex]
	write.table(genoRaw, file.path(resPath, "genotype", paste0("genotype_", chr, ".txt")), sep = "\t", quote = FALSE)
	write.table(genoMap, file.path(resPath, "genotype", paste0("genotype_", chr, "_Info.txt")), sep = "\t", quote = FALSE)
	
	#match meth data
	pIndex<-which(probeAnnot$chrm == chr)
	meth<-celltypeNormbeta[pIndex,]
	write.table(meth, file.path(resPath, "methylation", paste0("methylation_", chr, ".txt")), sep = "\t", quote = FALSE)
	write.table(probeAnnot[pIndex,c("probeID", "chrm", "start", "end")], file.path(resPath, "methylation", paste0("methylation_", chr, "_Info.txt")), sep = "\t", quote = FALSE)	
}


#----------------------------------------------------------------------#
# PROCESS COVARIATES
#----------------------------------------------------------------------#

genoPCs<-read.table(file.path(genoDir, "3_imputed","ImputationOutput", pop, "hg38", "MergedBatches_rsq0.3_QCd.pca.eigenvec"), stringsAsFactors = FALSE)
iIndex<-match(QCmetrics$Genotype.IID, idMap$V4[match(genoPCs$V1, idMap$V3)])
genoPCs<-genoPCs[iIndex,]

##PC plot
pdf(file.path(resPath, "PCPlotIncludedSamples.pdf"), width = 10, height = 10)
par(mfrow = c(2,2))
par(mar = c(4,4,0.5,0.5))
plot(genoPCs[,3], genoPCs[,4], pch = 16, xlab = "PC1", ylab = "PC2")
plot(genoPCs[,3], genoPCs[,5], pch = 16, xlab = "PC1", ylab = "PC3")
plot(genoPCs[,3], genoPCs[,6], pch = 16, xlab = "PC1", ylab = "PC4")
plot(genoPCs[,3], genoPCs[,7], pch = 16, xlab = "PC1", ylab = "PC5")
dev.off()

# combine continuous covariates
contCov<-cbind(QCmetrics$CCDNAmAge, genoPCs[,3:7])

# create dummy variable for sex
sex_m<-ifelse(QCmetrics$Sex == "M", 1, 0)

# create dummy variables for cell types
ct_neun<-ifelse(QCmetrics$Cell.type == "NeuN+", 1, 0)
ct_sox10<-ifelse(QCmetrics$Cell.type == "Sox10+", 1, 0)

## for interaction need covariate files with different cell types at the end
write.table(t(cbind(contCov, sex_m, ct_neun, ct_sox10)), file.path(resPath, "covariate", "covariates_sox10.txt"), sep = "\t", quote = FALSE)	
write.table(t(cbind(contCov, sex_m, ct_neun, ct_sox10)), file.path(resPath, "covariate", "covariates_neun.txt"), sep = "\t", quote = FALSE)	
