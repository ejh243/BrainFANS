## Written by Eilis
## summarize CNV QC metrics across samples
## provide path to sample qc output (from PennCNV) on command line
## determine thresholds

args<-commandArgs(trailingOnly = TRUE)
fileName<-args[1]
superPop<-args[2]
folder<-dirname(fileName)

## outlier defined as more than 3 SD from mean
nSD<-3


dat<-read.table(fileName, header = TRUE, stringsAsFactors = FALSE)

## use absolute WF
dat$absWF<-abs(dat$WF)

colsToPlot<-c("LRR_SD","BAF_drift","absWF","NumCNV")

sys.variables<-NULL
pdf(paste(folder, "HistogramCNVMetrics.pdf", sep = "/"))
par(mfrow = c(2,2))
par(mar = c(4,4,1,1))
for(colName in colsToPlot){

	hist(dat[,colName], xlab = colName, breaks = 20, main = "")
	mu<-mean(dat[,colName])
	sigma<-sd(dat[,colName])
	abline(v=mu, col = "red")

	abline(v=mu+nSD*sigma, col = "red", lty = 2)
	abline(v=mu-nSD*sigma, col = "red", lty = 2)
	# save thresholds for use with PennCNV filter functions
	sys.variables<-c(sys.variables, mu+nSD*sigma)
}

invisible(dev.off())

sys.variables<-signif(sys.variables, 3)

cat(sys.variables) # print so can be saved as bash variables


## create exclusion list
#exclude<-NULL
#for(colName in colsToPlot){
	#mu<-mean(dat[,colName])
	#sigma<-sd(dat[,colName])
	#exclude<-c(exclude, which(dat[,colName] > mu+nSD*sigma))
	#exclude<-c(exclude, which(dat[,colName] < mu-nSD*sigma))
#}
#exclude<-unique(exclude)
#print(paste(length(exclude), "samples excluded for outlier summary statistics"))

#gsub("PennCNVInput/", "", dat[exclude,1])
	
## compare metrics
pdf(paste(folder, "ScatterplotCNVMetrics.pdf", sep = "/"), height = 8, width = 12)
par(mar = c(4,4,1,1))
par(mfrow = c(2,3))
for(i in 1:3){
	for(j in c((i+1):4)){
		plot(dat[,colsToPlot[i]], dat[,colsToPlot[j]], pch = 16, xlab = colsToPlot[i], ylab = colsToPlot[j])
	}
}
invisible(dev.off())

## check by ethnicity
predPop<-read.csv(superPop)
predPop<-predPop[match(gsub("PennCNVInput/", "", dat$File), predPop$V1),]
pdf(paste(folder, "BoxplotCNVMetricsByPredPop.pdf", sep = "/"))
par(mfrow = c(2,2))
par(mar = c(4,4,1,1))
for(colName in colsToPlot){

boxplot(dat[,colName] ~ predPop$predPop, ylab = colName, main = "", col = rainbow(5))
}
invisible(dev.off())


	