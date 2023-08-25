if (! exists("project")){
	project <- ''
	dataType <- ''
} else if (! exists("dataType")){
	dataType <- ''	
}

dataDir<- paste0("/lustre/projects/Research_Project-MRC190311/", project)
metaDir<-paste0(dataDir, "/0_metadata")


chromDir<- paste0("/lustre/projects/Research_Project-MRC190311/integrative/chromHMM/", intproject)
completeDir<-paste0(chromDir,"/", dataType)
inDir<-paste0(completeDir, "/1_input")
tmpDir<-paste0(completeDir, "/2_binarised")
outDir<-paste0(completeDir, "/3_output")


binariseDir<-paste0(chromDir, "/1_binarised")
mergeDir<-paste0(chromDir, "/2_mergeBinarised")
modelDir<- paste0(chromDir, "/3_model")
qcDir<-paste0(modelDir, "/QCOutput" )
logDir<-"/lustre/home/jms260/BrainFANS/integrative/chromHMM/logFiles/jms260"
fastQCDir<-paste0(dataDir, "/1_raw/fastqc")
alignedDir<-paste0(dataDir, "/3_aligned")
