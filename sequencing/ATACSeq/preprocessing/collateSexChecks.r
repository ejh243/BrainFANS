#
args<-commandArgs(trailingOnly=TRUE)
dataPath<-args[1]

library(ggplot2)
library(ggpubr)
library(tidyverse)


## reads in peaks annotated to XIST and FIRRE
sampleSheet<-read.csv(file.path(dataPath, "0_metadata", "sampleSheet.csv"), stringsAsFactors = FALSE)

xPeakFile<-file.path(dataPath, "5_countPeaks", "MACS","ShiftedTagAlign", "sexChr", "chrX.peakcounts.txt")
yPeakFile<-file.path(dataPath, "5_countPeaks", "MACS","ShiftedTagAlign","sexChr", "chrY.peakcounts.txt")
xPeaks<-read.table(xPeakFile, stringsAsFactors = FALSE)
yPeaks<-read.table(yPeakFile, stringsAsFactors = FALSE)


xPeaks$V18<-basename(xPeaks$V18)
yPeaks$V10<-basename(yPeaks$V10)
xPeaksWide = xPeaks %>% 
	select(V4, V17, V18) %>% 
  spread(V4, V17)

yPeaksWide = yPeaks %>% 
	select(V4, V10, V11) %>% 
  spread(V4, V11)
  
xPeaksWide$YTot<-rowSums(yPeaksWide[,-1])
xPeaksWide$sex<-sampleSheet$gender[match(gsub(".chrX.tn5.tagAlign.gz", "", xPeaksWide$V18), sampleSheet$sampleID)]
 

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

figa<-ggarrange(a1, a2,  
          ncol = 2, nrow = 1)

pdf(file.path(dataPath, "3_aligned", "QCOutput", "sexPredictions.pdf"), width = 10, height = 5)
print(figa)
dev.off()
		  
## predict sex

sexPredict<-data.frame("sampleID" = gsub(".chrX.tn5.tagAlign.gz", "", xPeaksWide$V18),"LabelledSex" = xPeaksWide$sex, "FIRREratio" = ratioF > thresF,"XISTratio" = ratioX > thresX)
## TRUE means male
sexPredict[sexPredict == TRUE]<-"M"
sexPredict[sexPredict == FALSE]<-"F"
sexPredict$concordantPredictions<-sexPredict$FIRREratio == sexPredict$XISTratio
sexPredict$concordantLabel<-xPeaksWide$sex == sexPredict$FIRREratio
sexPredict$concordantLabel[sexPredict$concordantPredictions != TRUE]<-FALSE

write.csv(sexPredict, file.path(dataPath, "3_aligned", "QCOutput", "sexPredictions.csv"))