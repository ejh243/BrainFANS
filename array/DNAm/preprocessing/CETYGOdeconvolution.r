##---------------------------------------------------------------------#
##
## Title: Run CETYGO cell type deconvolution
##
## Purpose of script: Use CETYGO brain ref panal to run CETYGO on the samples
##
## Author: Emma Walker
##
## Date Created: Feb 2024
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# project folder is provided on command line
# More info on the CETYGO package can be found here: 
#   https://github.com/ejh243/CETYGO/wiki/Deconvolution-of-brain-cell-types


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

gdsFile <-paste0(dataDir, "/2_gds/raw.gds")
configFile <- paste0(dataDir, "/config.r")

gdsObj<-ifelse(file.exists(gdsFile), TRUE, ifelse(file.exists(msetFile), FALSE, NA))

source(configFile)

arrayType <- toupper(arrayType)


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(CETYGO)
library(gridExtra)


if(is.na(gdsObj)){
	message("No data to load")
} else{
	if(gdsObj){
		message("Loading data from gds object")
		library(bigmelon)
	} else {
		if(!gdsObj){
			message("Loading data from mset object")	
			library(wateRmelon)
		}
}
}



#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

setwd(dataDir)

# load sample sheet
sampleSheet<-read.csv("0_metadata/sampleSheet.csv", na.strings = c("", "NA"), stringsAsFactors = FALSE)
# if no column Basename, creates from columns Chip.ID and Chip.Location
if(!"Basename" %in% colnames(sampleSheet)){
	sampleSheet$Basename<-paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
}

if(gdsObj){

	gfile<-openfn.gds(gdsFile, readonly = FALSE, allow.fork = TRUE)
	# ensure sample sheet is in same order as data
	sampleSheet<-sampleSheet[match(colnames(gfile), sampleSheet$Basename),]
	# extract a few useful matrices
	if(arrayType == "V2"){
	  rawbetas<-epicv2clean(betas(gfile)[])
} else {
	rawbetas<-betas(gfile)[,]
}
    closefn.gds(gfile)
}else {
	if(!gdsObj){
		load(msetFile)
		# ensure sample sheet is in same order as data
		sampleSheet<-sampleSheet[match(colnames(msetEPIC), sampleSheet$Basename),]
		# extract a few useful matrices
		rawbetas<-betas(msetEPIC)
}
}


  


#----------------------------------------------------------------------#
# DEFINE CETYGO FUNCTION
#----------------------------------------------------------------------#

adultBrainCETYGO <- function(betas, cellType){
  
  predPropAll<-list()
  CETYGOplots <- list()
  propPlots <- list()
  counter = 1
  maxCETYGO <- 0
  minCETYGO <- 1
  
  for(method in names(modelBrainCoef)){
    for(j in 1:length(modelBrainCoef[[method]])){
      if(!is.null(modelBrainCoef[[method]][[j]])){
        predPropAll[[method]][[j]]<-projectCellTypeWithError(betas, modelBrainCoef[[method]][[j]])
      }
    } 
  }
  
  
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[1]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        if (max(plotdf$CETYGO) > maxCETYGO){
          maxCETYGO <- max(plotdf$CETYGO)}
        
        if (min(plotdf$CETYGO) < minCETYGO){
          minCETYGO <- min(plotdf$CETYGO)}
      }
    }
  }
  
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[1]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        # CETYGO score boxplot
        pCETYGO <- ggplot(plotdf, aes(factor(0), CETYGO))+
          geom_boxplot()+
          coord_cartesian(ylim=c(minCETYGO - 0.01, maxCETYGO + 0.01))+
          ggtitle(paste0(cellType, "_", names(predPropAll[i]), "-", j))+
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.title = element_text(size=6))
        
        CETYGOplots[[counter]] <- pCETYGO
        
        # proportions boxplot
        plotdf <-as.data.frame(plotdf[,1:(ncol(plotdf)-2)])
        plotdf$Basename <- rownames(plotdf)
        plotdf <- reshape2::melt(plotdf, id.vars = "Basename")
        colnames(plotdf) <- c("Basename", "CellType", "Proportion")
        
        p <- ggplot(plotdf, aes(x=CellType, y = Proportion)) +
          geom_boxplot()+
          ggtitle(paste0(cellType, "_", names(predPropAll[i]), "-", j))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        propPlots[[counter]] <- p
        
        counter = counter + 1
      }
    }
  }
  
  CETYGOplotFile <- paste0(dataDir, "/2_gds/QCmetrics/CETYGO/CETYGOscorePlots_", cellType, ".pdf")
  pdf(CETYGOplotFile)
  do.call(grid.arrange, c(CETYGOplots, ncol = 4))
  dev.off()
  
  propPlotsFile <- paste0(dataDir, "/2_gds/QCmetrics/CETYGO/CETYGOcellPropsPlots_", cellType, ".pdf")
  pdf(propPlotsFile, height = 20, width = 20)
  do.call(grid.arrange, c(propPlots, ncol = 4))
  dev.off()
  
  predPropOutFile <- paste0(dataDir, "/2_gds/QCmetrics/CETYGO/CETYGOpredictedProportions_", cellType, ".rdat")
  save(predPropAll, file = predPropOutFile)
  
}


#----------------------------------------------------------------------#
# RUN CETYGO FUNCTION ON EACH CELL TYPE
#----------------------------------------------------------------------#

for(cell in sampleSheet$Cell_Type){
    cellSampleSheet <- sampleSheet[which(sampleSheet$Cell_Type == cell),]
    cellBetas <- rawbetas[, colnames(rawbetas) %in% cellSampleSheet$Basename]
    adultBrainCETYGO(cellBetas, cell)
}
