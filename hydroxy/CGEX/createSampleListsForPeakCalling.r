## create lists of samples for peak calling by sample type
args<-commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("Config file missing on command line", call.=FALSE)
} else {
  print(paste("Using the following config file:", args[1]))
}

source(args[1]) ## contains paths to files

setwd(dataDir)
pheno<-read.csv(sampleSheet, stringsAsFactors = FALSE, row.names = 1)

## reformat filenames to just basenames
pheno$bamReads<-gsub("alignedData/", "", pheno$bamReads)
pheno$controlReads<-gsub("alignedData/", "", pheno$controlReads)

print("Sample sheet loaded")

## by brain region
for(each in unique(pheno$Tissue)){
	index<-which(pheno$Tissue == each)
	input<-pheno$bamReads[index]
	controls<-pheno$controlReads[index]
	write.table(input, paste0("sampleLists/PeakCallingInputFiles", each, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
	write.table(controls, paste0("sampleLists/PeakCallingControlFiles", each, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

print("Brain Regions done")

## by cell type
for(each in unique(pheno$Fraction)){
	index<-which(pheno$Fraction == each)
	input<-pheno$bamReads[index]
	controls<-pheno$controlReads[index]
	write.table(input, paste0("sampleLists/PeakCallingInputFiles", each, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
	write.table(controls, paste0("sampleLists/PeakCallingControlFiles", each, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

print("Fractions done")

## by cell type x fraction
for(each in unique(pheno$Tissue)){
	for(item in unique(pheno$Fraction)){
		index<-which(pheno$Tissue == each & pheno$Fraction == item)
		input<-pheno$bamReads[index]
		controls<-pheno$controlReads[index]
		write.table(input, paste0("sampleLists/PeakCallingInputFiles", each,"_",item, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
		write.table(controls, paste0("sampleLists/PeakCallingControlFiles", each,"_",item, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
	}
}

print("Brain Regions x Fractions done")


