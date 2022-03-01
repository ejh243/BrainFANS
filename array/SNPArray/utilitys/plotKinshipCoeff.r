## plot kinship coefficients

args<-commandArgs(TRUE)

datFileName <- args[1]

if(dir.exists(datFileName)){
	allFiles <- list.files(datFileName, pattern = "king.kin0")
	dat <- NULL
	for(file in allFiles){
		dat.tmp<-read.table(paste0(datFileName, file), header = TRUE, stringsAsFactors = FALSE)
		dat<-rbind(dat, dat.tmp)
	}
	population<-"All"
	outFolder<-datFileName
} else {
	population <- lapply(strsplit(gsub("_QCd_king.kin0", "", basename(datFileName)), "_"), tail , n = 1)[[1]]
	dat<-read.table(datFileName, header = TRUE, stringsAsFactors = FALSE)
	outFolder<-dirname(datFileName)
}

pdf(paste0(outFolder, "/ScatterplotKinshipCoefficients", population, ".pdf"))
plot(dat$IBS0, dat$Kinship, pch = 16,
  col = "blue", main = "KING-Robust (Default)",
  xlab="Proportion of Zero IBS", ylab = "Estimated Kinship Coefficient", ylim = c(0,0.4))
  abline(h = 0.3536, col = "black", lty = 3)
abline(h = 0.1768, col = "black", lty = 3)
abline(h = 0.0884, col = "black", lty = 3)
abline(h = 0.0442, col = "black", lty = 3)
abline(h = 0.0221, col = "black", lty = 3)
  dev.off()
  
## look for samples with higher than average IBD statistics
sampleIDs<-unique(c(dat$FID1, dat$FID2))
mean.kin<-rep(NA, length(sampleIDs))
names(mean.kin)<-sampleIDs
for(each in sampleIDs){
	values<-c(dat$Kinship[which(dat$FID1 == each | dat$FID2 == each)])
	mean.kin[each]<-mean(values, na.rm = TRUE)
}
mean.ibs<-rep(NA, length(sampleIDs))
names(mean.ibs)<-sampleIDs
for(each in sampleIDs){
	values<-c(dat$IBS0[which(dat$FID1 == each | dat$FID2 == each)])
	mean.ibs[each]<-mean(values, na.rm = TRUE)
}


pdf(paste0(outFolder, "/HistMeanKinshipCoefficients", population, ".pdf"))
par(mfrow = c(1,2))
hist(mean.kin, breaks = 20, xlab = "Mean kinship coefficient", main = "")
mu<-mean(mean.kin)
sigma<-sd(mean.kin)
abline(v = mu, col = "red", lty = 2)
for(i in 1:3){
	abline(v = mu+i*sigma, col = "red")
	abline(v = mu-i*sigma, col = "red")
}
hist(mean.ibs, breaks = 20, xlab = "Mean IBS coefficient", main = "")
mu<-mean(mean.ibs)
sigma<-sd(mean.ibs)
abline(v = mu, col = "red", lty = 2)
for(i in 1:3){
	abline(v = mu+i*sigma, col = "red")
	abline(v = mu-i*sigma, col = "red")
}

## look for joint outliers

par(mfrow = c(1,1))
plot(mean.ibs, mean.kin, pch = 16)
mu<-mean(mean.ibs)
sigma<-sd(mean.ibs)
abline(v = mu, col = "red", lty = 2)
for(i in 1:3){
	abline(v = mu+i*sigma, col = "red")
	abline(v = mu-i*sigma, col = "red")
}
mu<-mean(mean.kin)
sigma<-sd(mean.kin)
abline(h = mu, col = "red", lty = 2)
for(i in 1:3){
	abline(h = mu+i*sigma, col = "red")
	abline(h = mu-i*sigma, col = "red")
}
dev.off()

#excessRelated<-names(which((mean.ibs > mean(mean.ibs)+i*sd(mean.ibs) | mean.ibs < mean(mean.ibs)-i*sd(mean.ibs))  & (mean.kin >  mean(mean.kin)+i*sd(mean.kin) | mean.kin <  mean(mean.kin)-i*sd(mean.kin))))

#dat[dat$FID1 %in% excessRelated | dat$FID2 %in% excessRelated,]