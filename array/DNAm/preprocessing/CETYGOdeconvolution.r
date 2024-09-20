##---------------------------------------------------------------------#
##
## Title: Run CETYGO cell type deconvolution
##
## Purpose of script: Use CETYGO brain ref panal or blood ref panal depending on ## tissue type to run CETYGO on the samples
##
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

source(configFile)

arrayType <- toupper(arrayType)
tissueType <- toupper(tissueType)


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(CETYGO)
library(gridExtra)
library(bigmelon)
library(ggplot2)
library(wateRmelon)


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
sampleSheet$Cell_Type <- as.factor(sampleSheet$Cell_Type)


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



#----------------------------------------------------------------------#
# DEFINE CETYGO FUNCTION FOR BRAIN TISSUE
#----------------------------------------------------------------------#

adultBrainCETYGO <- function(betas, cellType){
  
  predPropAll<-list() # store output in list
  CETYGOplots <- list() # store CETYGO score plots in list
  propPlots <- list() # store proportion plots in list
  counter = 1 # counter for storing list elements from for loop
  maxCETYGO <- 0 # vector to hold max CETYGO score from each model/ref
  minCETYGO <- 1 # vector to hold min CETYGO score from each model/ref
  
  # run CETYGO using both ANOVA and IDOL methods and all available cell ref panals 
  for(method in names(modelBrainCoef)){
    for(j in 1:length(modelBrainCoef[[method]])){
      if(!is.null(modelBrainCoef[[method]][[j]])){
        predPropAll[[method]][[j]]<-projectCellTypeWithError(betas, modelBrainCoef[[method]][[j]])
      }
    } 
  }
  
  # Find max and min CETYGO scores across all methods/ref panals to calibrate axis limits
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[i]]
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
  
  # create boxplots of CETYGO scores and proportions
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[i]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        # CETYGO score boxplot
        pCETYGO <- ggplot(plotdf, aes(factor(0), CETYGO))+
          geom_boxplot()+
          coord_cartesian(ylim=c(minCETYGO - 0.01, maxCETYGO + 0.01))+
          ggtitle(paste0(cellType, "_", names(predPropAll[i]), "-", j))+
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.title = element_text(size=6))+
          geom_hline(yintercept=0.07, linetype="dashed", color = "red")
        
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
  
  # arrange plots into to single object and save all output
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
# DEFINE CETYGO FUNCTION FOR BLOOD TISSUE
#----------------------------------------------------------------------#


adultBloodCETYGO <- function(betas){
  
  # run CETYGO using blood ref panal
  rowIndex<-rownames(betas)[rownames(betas) %in% rownames(modelBloodCoef)]
  predProp<-as.data.frame(projectCellTypeWithError(betas, modelBloodCoef[rowIndex,]))


  # boxPlot CETYGO score
  pCETYGO <- ggplot(predProp, aes(factor(0), CETYGO))+
    geom_boxplot()+
  geom_hline(yintercept=0.07, linetype="dashed", color = "red")+
    xlab("")
  
  #save plot
  CETYGOplotFile <- paste0(dataDir, "/2_gds/QCmetrics/CETYGO/CETYGOscorePlot.pdf")
  pdf(CETYGOplotFile)
  print(pCETYGO)
  dev.off()

  
  # box plot proportions
  plotdf <-as.data.frame(predProp[,1:(ncol(predProp)-2)])
  plotdf$Basename <- rownames(plotdf)
  plotdf <- reshape2::melt(plotdf, id.vars = "Basename")
  colnames(plotdf) <- c("Basename", "CellType", "Proportion")

  p <- ggplot(plotdf, aes(x=CellType, y=Proportion)) +
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # save plot
  propPlotsFile <- paste0(dataDir, "/2_gds/QCmetrics/CETYGO/CETYGOcellPropsPlot.pdf")
  pdf(propPlotsFile)
  print(p)
  dev.off()
  
  # save output file
  predPropOutFile <- paste0(dataDir, "/2_gds/QCmetrics/CETYGO/CETYGOpredictedProportions.rdat")
  save(predProp, file = predPropOutFile)

}



#----------------------------------------------------------------------#
# RUN CETYGO FUNCTION 
#----------------------------------------------------------------------#

# for sorted Brain tissue run on each cell type individually
if(tissueType == "BRAIN" && ctCheck) {

  for(cell in levels(sampleSheet$Cell_Type)) {
    cellSampleSheet <- sampleSheet[which(sampleSheet$Cell_Type == cell), ]
    cellBetas <- rawbetas[, colnames(rawbetas) %in% cellSampleSheet$Basename]
    adultBrainCETYGO(cellBetas, cell)
  }
} else {
  if(tissueType == "BLOOD") {
    adultBloodCETYGO(rawbetas)
  }
  if(tissueType == "BRAIN") {
    adultBrainCETYGO(rawbetas, "bulk")
  }
}
