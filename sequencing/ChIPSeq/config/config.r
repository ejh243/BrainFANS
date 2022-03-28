## Sequencing data filepaths
dir<- paste0("/gpfs/mrc0/projects/Research_Project-MRC190311/ChIPSeq/", project)

metaDir<- paste0(dir, "/0_metadata")
dataDir<-paste0(dir, "/1_raw")
fastQCDir<-paste0(dataDir, "/fastqc")
trimDir<-paste0(dir, "/2_trimmed") 
alignedDir<-paste0(dir, "/3_aligned")
peakDir<- paste0(dir, "/4_calledPeaks")
qcDir<-paste0(alignedDir, "/QCOutput" )
sampleSheet<-paste0(metaDir, "/sampleSheet.csv")
