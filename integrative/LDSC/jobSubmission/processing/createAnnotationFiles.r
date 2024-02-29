## You need BEDTools and BEDOPs available on the command line for this script
## to work properly.

configuration_file_path <- commandArgs(trailingOnly = TRUE)

source(configuration_file_path)

library(data.table)
library(bedr)
library(dplyr)


read_peaks <- function(file_paths, peaks_names) {
	if (length(file_paths) != length(peaks_names)) {
		stop("number of file paths does not match number of peak names given")
	}
	sorted_peaks <- list()

	for (peak in seq_along(list_of_peaks)) {
		# Sorts by chromosome and start position in case the peak calls are
		# not already in this format. Sorting speeds up bedr functions.
		peaks <- read.table(file_paths[peak], stringsAsFactors = FALSE)
		peaks <- peaks[order(peaks[, 1], peaks[, 2], )]

		regions <- bed2index(peaks[, 1:3])
		sorted_peaks[[peak]] <- regions
	}
	return(sorted_peaks)
}

sorted_peaks <- read_peaks(peaks_file_paths, peaks_names)

## This bedr function is expecting one bed file, so we bind all of the peak
## files into one long bed file. The sort is for faster bedr processing
## for in.region.
## Commented out as not currently in use.
# all_regions <- bedr.merge.region(bedr.sort.region(bind_rows(sorted_peaks)))

for(chr in 1:22) {
	annot <-read.table(gzfile(paste0(baseline_file_prefix, chr, ".annot.gz")),
										header = TRUE)

	snps<-fread(paste0(bim_file_prefix, chr, ".bim"))

	## This is in place to remove any annotation columns that are irrelevant
	## or non-binary. Binary columns are required later for determining if a 
	## snp resides in any of the baseline regulatory regions.
	annot <- annot %>% select_if(function(x) all(is.integer(x) | x %in% c(0,1)))
	
	## We extract the base pair position twice to get a start
	## and end read for each SNP (SNPs are one base pair always, so we only need
	## the start BP and add one on to this value to get the end).
	snp_bed<-annot[,c("CHR", "BP", "BP")]
	snp_bed<-snp_bed[order(snp_bed$CHR, snp_bed$BP),]
	colnames(snp_bed)<-c("chr", "start", "end")
	snp_bed$end<-snp_bed$end+1
	snp_bed$chr<-paste("chr", snp_bed$chr, sep = "")
	class(snp_bed$end)<-"integer"
	snp_regions<-bed2index(as.data.frame(snp_bed))
	
	## classify SNPs in the called peaks
	snp_classifications <-
	  lapply(sorted_peaks, function(peaks) in.region(snp_regions, peaks))
	
	## This is unused in the output and so is commented out. Adding this in
	## is liable to introduce colinearity between annoation columns
	# overlapany <- in.region(snpRegions, all_regions)

	mergeOverlaps<-bind_cols(snp_classifications)

	## Creates a column for snps that are in any feature given in the baseline
	## annotation files.
	## We exclude columns one to three as these columns are always positive
	## (chromosome, position and base).  
	any_regulatory_feature <- as.numeric(rowSums(annot[,-c(1:3)]) > 0)

	output_annotation <-
  data.frame(snps$V1, snps$V4, snps$V2, snps$V3,
             mergeOverlaps, any_regulatory_feature, annot$base)

	colnames(output_annotation) <-
	  c("CHR", "BP", "SNP", "CM", peaks_names, "AnyRegulatoryFeature", "AllSNPs")

	## This is to ensure R does not add quotes around integers, which would 
	## result in errors once parsed into ldsc
	output_annotation[peaks_names] <-
	  lapply(output_annotation[peaks_names], as.numeric)

	annotation_file_name <-
	  paste0(ld_annotation_dir, "/", ld_annotation_prefix, ".", chr, ".annot.gz")

	write.table(output_annotation,
            annotation_file_name,
            quote = FALSE,
            row.names = FALSE)
}