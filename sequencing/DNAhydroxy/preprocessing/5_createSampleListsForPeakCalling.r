## create lists of samples for peak calling by sample type
args<-commandArgs(trailingOnly = TRUE)
dataDir<-args[1] ## contains paths to files

alignedPath<-file.path(dataDir, "3_aligned")


setwd(dataDir)
pheno<-read.csv(file.path("0_metadata", "sampleSheet.csv"), stringsAsFactors = FALSE, row.names = 1)

## reformat filenames to just basenames
pheno$bamReads<-gsub("3_aligned/", "", pheno$bamReads)
pheno$controlReads<-gsub("3_aligned/", "", pheno$controlReads)

print("Sample sheet loaded")

## file for all samples
write.table(pheno$bamReads, file.path(dataDir, "0_metadata", "sampleLists", "PeakCallingInputFilesALL.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(pheno$controlReads, file.path(dataDir, "0_metadata", "sampleLists", "PeakCallingControlFilesALL.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

## by brain region
for(each in unique(pheno$tissue)){
	index<-which(pheno$tissue == each)
	input<-pheno$bamReads[index]
	controls<-pheno$controlReads[index]
	write.table(input, file.path(dataDir, "0_metadata", "sampleLists", paste0("PeakCallingInputFiles", each, ".txt")), quote = FALSE, col.names = FALSE, row.names = FALSE)
	write.table(controls, file.path(dataDir, "0_metadata", "sampleLists", paste0("PeakCallingControlFiles", each, ".txt")), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

print("Brain Regions done")

## by cell type
for(each in unique(pheno$F)){
	index<-which(pheno$fraction == each)
	input<-pheno$bamReads[index]
	controls<-pheno$controlReads[index]
	write.table(input, file.path(dataDir, "0_metadata", "sampleLists", paste0("PeakCallingInputFiles", each, ".txt")), quote = FALSE, col.names = FALSE, row.names = FALSE)
	write.table(controls, file.path(dataDir, "0_metadata", "sampleLists", paste0("PeakCallingControlFiles", each, ".txt")), quote = FALSE, col.names = FALSE, row.names = FALSE)

}

print("fractions done")

## by cell type x fraction
for(each in unique(pheno$tissue)){
	for(item in unique(pheno$fraction)){
		index<-which(pheno$tissue == each & pheno$fraction == item)
		input<-pheno$bamReads[index]
		controls<-pheno$controlReads[index]
		write.table(input, file.path(dataDir, "0_metadata", "sampleLists", paste0("PeakCallingInputFiles", each,"_",item, ".txt")), quote = FALSE, col.names = FALSE, row.names = FALSE)
		write.table(controls, file.path(dataDir, "0_metadata", "sampleLists", paste0("PeakCallingControlFiles", each,"_",item, ".txt")), quote = FALSE, col.names = FALSE, row.names = FALSE)
	}
}

print("Brain Regions x fractions done")


