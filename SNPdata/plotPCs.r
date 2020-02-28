## saves output to current folder

args<-commandArgs(TRUE)

filename<-args[1]
nSD<-as.numeric(args[2])

impute<-read.table(filename)

pcMu<-colMeans(impute[,-c(1:2)])
pcSigma<-apply(impute[,-c(1:2)], 2, sd)

pdf("ScatterplotPCs.pdf", width = 12, height = 4)
par(mfrow = c(1,3))
par(mar = c(4,4,1,1))
plot(impute[,3],impute[,4], pch= 16, xlab = "PC1", ylab = "PC2")
abline(v = pcMu[1]+nSD*pcSigma[1], lty = 2)
abline(v = pcMu[1]-nSD*pcSigma[1], lty = 2)
abline(h = pcMu[2]+nSD*pcSigma[2], lty = 2)
abline(h = pcMu[2]-nSD*pcSigma[2], lty = 2)

plot(impute[,5],impute[,6], pch= 16, xlab = "PC3", ylab = "PC4")
abline(v = pcMu[3]+nSD*pcSigma[3], lty = 2)
abline(v = pcMu[3]-nSD*pcSigma[3], lty = 2)
abline(h = pcMu[4]+nSD*pcSigma[4], lty = 2)
abline(h = pcMu[4]-nSD*pcSigma[4], lty = 2)

plot(impute[,7],impute[,8], pch= 16, xlab = "PC5", ylab = "PC6")
abline(v = pcMu[5]+nSD*pcSigma[5], lty = 2)
abline(v = pcMu[5]-nSD*pcSigma[5], lty = 2)
abline(h = pcMu[6]+nSD*pcSigma[6], lty = 2)
abline(h = pcMu[6]-nSD*pcSigma[6], lty = 2)
dev.off()


## write list of outliers
outliers<-NULL
for(i in 3:22){
  index<-which(impute[,i] > pcMu[i-2]+nSD*pcSigma[i-2] | impute[,i] < pcMu[i-2]-nSD*pcSigma[i-2])
  if(length(index) > 0){
   outliers<-rbind(outliers, cbind(impute[index,1:2], paste("Outlier PC", (i-2), sep = "")))
  }
}

write.table(outliers, paste("OutliersFromPC_",nSD,"SDfromMean.txt", sep = ""), quote = FALSE, sep = " ")