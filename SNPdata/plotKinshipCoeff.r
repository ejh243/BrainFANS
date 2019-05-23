## plot kinship coefficients

args<-commandArgs(TRUE)

setwd(args[1])

dat<-read.table("king.kin0", header = TRUE, stringsAsFactors = FALSE)
outPred<-read.csv("PredictedPopulations.csv", stringsAsFactors = FALSE)


pdf("ScatterplotKinshipCoefficients.pdf")
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

outPred<-outPred[match(names(mean.ibs), outPred$V1),]

pdf("HistMeanKinshipCoefficients.pdf")
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
## colour by predicted super population

par(mfrow = c(1,1))
plot(mean.ibs, mean.kin, pch = 16, col = rainbow(5)[as.factor(outPred$predPop)])
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
