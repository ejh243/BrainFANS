##---------------------------------------------------------------------#
##
## Title: cellMarkBamFile
##
## Purpose of script: to generate a txt file in the following format (tsv), named cellMarkFileTable.txt
##                                    cell1 	mark1 	case1.bam 	control1.bam
##                                    cell2	  mark2 	case2.bam 	control2.bam
##
## Author: Jessica Shields
##
## Date Created: 2022-07-05
##
##---------------------------------------------------------------------#
##
## clear the R environment
rm(list=ls()) 

## set working directory
setwd("/lustre/projects/Research_Project-MRC190311/scripts/sequencing")

## load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args[1]<-"DNAhydroxy/MRC"
  #args[1]<-'WGBS/rizzardi'
  args[1]<-'ATACSeq/rizzardi'
  args[2]<-'prefrontal cortex|PFC'
} 

project<-args[1]
tissue<-args[2]
source("BSSeq/config/config.r")
#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
sampleSheet<-read.csv(sampleSheet)

if (file.exists(paste0(metaDir, '/stage1Samples.txt'))){
  sampleid<-read.table(paste0(metaDir, '/stage1Samples.txt'))[,1]
} else {
  print('Using unfiltered samples.txt')
  sampleid<-read.table(paste0(metaDir, '/samples.txt'))[,1]
}

#----------------------------------------------------------------------#
# GENERATE CELLMARKBAMFILE
#----------------------------------------------------------------------#
sampleSheet<-sampleSheet[which(sampleSheet$sampleID %in% sampleid),]

## filter sampleSheet by selected tissue
if (length(sampleSheet$tissue)!=0){
  sampleSheet<-sampleSheet[grep(tissue, sampleSheet$tissue), ]
}
sampleSheet


# create reference dictionary to rename fractions
dic<-data.frame( 
    c('N+', 'glu|gaba|neun|neun pos'),
    c('N-', 'olig|sox10|neun neg'),
    c('T', 'total|bulk'), 
    c('DN', 'neg')
    )
# rename fractions to N-, N- and S
cell <- tolower(sampleSheet$fraction)
for (x in 1:length(colnames(dic))){
  cell<-replace(cell, 
                grepl(dic[2,x], cell), 
                as.character(dic[1,x]))
}


## get assay target
mark<-tolower(sampleSheet$target)


##if files have bed files use them (atac) otherwise use bam
reads<-list.files(alignedDir, pattern = "gz$")
if (length(reads) != 0){
  index<-match(sampleSheet$sampleID, gsub('.tn5.tagAlign.gz', "", reads))
} else {
  ## get bam files (agnostic of data type) 
  reads<-list.files(alignedDir, pattern = ".\\wdup.+bam$")
  # endings differ so use sampleID to check length of file suffix
  bamString<-substr(reads[grep(paste0(sampleSheet$sampleID[1],'.+'), reads)],
                    nchar(levels(sampleSheet$sampleID)[1])+1, nchar(reads[1])+1)
  index<-match(sampleSheet$sampleID, gsub(bamString, "", reads))
}
reads<-reads[index]


## if files include a control, get list of control bams
if (length(sampleSheet$controlID)!=0){
  controlReads<-list.files(alignedDir, pattern = ".\\wdup.+bam$")
  index<-match(sampleSheet$controlID, gsub(bamString, "", controlReads))
  controlReads<-controlReads[index]
  
  cellMarkFile<-cbind(cell, mark, reads, controlReads)
} else {
  cellMarkFile<-cbind(cell, mark, reads)
}

# save to metaDir
write.table(cellMarkFile, paste0(metaDir, '/cellMarkFileTable.txt'), 
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
