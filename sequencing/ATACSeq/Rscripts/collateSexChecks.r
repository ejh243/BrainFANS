## ==========================================================================================##
##                      ATAC-seq pipeline STEP 5.3: Collate sex checks                       ##
## ==========================================================================================##
## EXECUTION: Rscript ./Rscripts/collateSexChecks.r <dataPath>                               ||
## - execute from scripts directory                                                          ||
##                                                                                           ||
## DESCRIPTION: This scripts collates results from calling peaks in sex chromosomes and      ||
## uses this to check the assigned sex of samples.                                           ||
##                                                                                           ||
## INPUTS:                                                                                   ||
## - <config file> : R config file for project                                               ||
##                                                                                           ||
## OUTPUTS:                                                                                  ||
## - sexPredictions.csv, sexPredictions.pdf                                                  ||
##                                                                                           ||
## REQUIRES:                                                                                 ||
## - R/4.2.1-foss-2022a, libraries: ggpubr, tidyverse,ggplot2                                ||
## ==========================================================================================##

args<-commandArgs(trailingOnly=TRUE)
configFile<-args[1]
source(configFile)

library(ggplot2)
library(ggpubr)
library(tidyverse)

## reads in peaks annotated to XIST and FIRRE
sampleSheet<-read.csv(file.path(dir, "/0_metadata", "sampleSheet.csv"), stringsAsFactors = FALSE,colClasses="character")
xPeakFile<-file.path(dir, "/5_countPeaks", "ShiftedTagAlign", "sexChr", "chrX.peakcounts.txt")
yPeakFile<-file.path(dir, "/5_countPeaks", "ShiftedTagAlign","sexChr", "chrY.peakcounts.txt")

xPeaks<-read.table(xPeakFile, stringsAsFactors = FALSE)
yPeaks<-read.table(yPeakFile, stringsAsFactors = FALSE,fill = TRUE)

colnames(xPeaks)[c(4,17,18)] <- c("sex-gene","counts","sampleID")
colnames(yPeaks)[c(4,10,11)] <- c("peak-name","counts","sampleID")

xPeaks$sampleID<-basename(xPeaks$"sampleID")
yPeaks$sampleID<-basename(yPeaks$"sampleID")

xPeaksWide = xPeaks %>% 
	select("sex-gene","counts","sampleID") %>% 
  spread("sex-gene","counts")

yPeaksWide = yPeaks %>% 
	select("peak-name","counts","sampleID") %>%
 spread("peak-name","counts")

xPeaksWide$sampleID<-gsub(".chrX.tn5.tagAlign.gz", "", xPeaksWide$sampleID, fixed = TRUE)
yPeaksWide$sampleID<-gsub(".chrY.tn5.tagAlign.gz", "", yPeaksWide$sampleID,fixed = TRUE)

xPeaksWide<-xPeaksWide[match(sampleSheet$sampleID,xPeaksWide$sampleID), ]
yPeaksWide<-yPeaksWide[match(sampleSheet$sampleID,yPeaksWide$sampleID),]

xPeaksWide$YTot<-rowSums(yPeaksWide[,-1])
xPeaksWide$sex<-sampleSheet$gender

## determine threshold as a line to bisect the two groups 
ratioF<-xPeaksWide$YTot/xPeaksWide$FIRRE
ratioF[!is.finite(ratioF)]<-NA
sexMeans<-aggregate(ratioF, by = list(xPeaksWide$sex), mean, na.rm = TRUE)[,2]
thresF<-min(sexMeans)+ 0.5*abs(diff(sexMeans))

ratioX<-xPeaksWide$YTot/xPeaksWide$XIST
ratioX[!is.finite(ratioX)]<-NA
sexMeans<-aggregate(ratioX, by = list(xPeaksWide$sex), mean, na.rm = TRUE)[,2]
thresX<-min(sexMeans)+ 0.5*abs(diff(sexMeans))

a1<-ggplot(xPeaksWide, aes(x=FIRRE, y=YTot, color=sex)) + 
  geom_point() + 
  labs(
    x = "Counts in FIRRE associated peak",
    y = "Counts across all Y chr peaks"
  ) + theme_bw() + geom_abline(intercept = 0, slope = thresF)
  
a2 <- ggplot(xPeaksWide, aes(x=XIST, y=YTot, color=sex)) + 
  geom_point() + 
  labs(
    x = "Counts in XIST associated peak",
    y = "Counts across all X chr peaks"
  ) + theme_bw() + geom_abline(intercept = 0, slope = thresX)

plots<-ggarrange(a1, a2,  
          ncol = 2, nrow = 1)

pdf(file.path(dir, "3_aligned", "QCOutput", "sexPredictions.pdf"), width = 10, height = 5)
print(plots)
dev.off()
		  
## predict sex
sexPredict<-data.frame("sampleID" = xPeaksWide$sampleID,"LabelledSex" = xPeaksWide$sex, "FIRREratio" = ratioF > thresF,"XISTratio" = ratioX > thresX)

## TRUE means male
sexPredict[sexPredict == TRUE]<-"M"
sexPredict[sexPredict == FALSE]<-"F"
sexPredict$concordantPredictions<-sexPredict$FIRREratio == sexPredict$XISTratio
sexPredict$concordantLabel<-xPeaksWide$sex == sexPredict$FIRREratio
sexPredict$concordantLabel[sexPredict$concordantPredictions != TRUE]<-FALSE

write.csv(sexPredict, file.path(dir, "/3_aligned", "QCOutput", "sexPredictions.csv"))