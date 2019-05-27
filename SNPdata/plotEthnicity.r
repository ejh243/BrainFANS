## script to plot PCA output and colour by 1000G ethnicity
args<-commandArgs(TRUE)

setwd(args[1])
pcas<-read.table("mergedw1000G.pca.eigenvec", stringsAsFactors = FALSE)
KGped<-read.table("../References/1000G/20130606_g1k.ped", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
popInfo<-read.table("../References/1000G/PopInfo.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t") ## table made from 1000G website

KGped<-KGped[match(pcas[,2], KGped[,2]),]
nPops<-length(table(KGped$Population))
popInfo<-popInfo[match(popInfo$Population.Code,levels(as.factor(KGped$Population))),]

KGped<-cbind(KGped,popInfo$Super.Population.Code[match(KGped$Population, popInfo$Population.Code)])
colnames(KGped)[ncol(KGped)]<-"SuperPopulation"
nSuperPops<-length(table(KGped$SuperPopulation))


ptType<-c(20,3)[as.factor(is.na(KGped[,1]))]
ptCol<-rainbow(nPops)[as.factor(KGped$Population)]
ptCol[is.na(ptCol)]<-"black"

pdf("PCAplotwith1KG.pdf", width = 10, height = 10)
par(mfrow = c(2,2))
par(mar = c(4,4,0.75,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nPops), levels(as.factor(KGped$Population)), cex = 1.5, ncol = 3)


## alternatively plot "super populations"


ptCol<-rainbow(nSuperPops)[as.factor(KGped$SuperPopulation)]
ptCol[is.na(ptCol)]<-"black"
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)
dev.off()


calcPopDist<-function(dat.pca, ref){
# dat.pca is test individual's data across pcs
# ref is matrix with 1 row per population and 1 column per pca
	popDist<-rep(NA, nrow(ref))
	names(popDist)<-rownames(ref)
	for(i in 1:nrow(ref)){
		diffs<-dat.pca - ref[i,]
		sqdiffs<-diffs^2
		popDist[i]<-sqrt(sum(sqdiffs))
	}
	return(popDist)
}

## for each super population calculate cluster medians

nMatches<-rep(NA,20) 
for(nPCs in 2:20){
	pop.medians<-apply(pcas[,-c(1:2)][,1:nPCs], 2,aggregate, by = list(KGped$SuperPopulation), median)
	pop.medians<-cbind.data.frame(pop.medians)
	rownames(pop.medians)<-pop.medians[,1]
	pop.medians<-pop.medians[,seq(2,nPCs*2,2)]

	## for each individual compare to each super population and find most similar


	popDistsAll<-apply(pcas[,-c(1:2)][,1:nPCs], 1, calcPopDist, pop.medians)
	popDistsAll<-t(popDistsAll)
	predPop<-colnames(popDistsAll)[apply(popDistsAll, 1, which.min)]

	compTrue<-table(predPop, KGped$SuperPopulation)
	nMatches[nPCs]<-sum(diag(compTrue))
}
pdf("SelectOptimalnPCsForPopulationPrediction.pdf")
plot(2:20,nMatches/sum(!is.na(KGped$SuperPopulation))*100, xlab = "nPCs", ylab = "Percentage Correct")
dev.off()

nPCs<-which.max(nMatches)
pop.medians<-apply(pcas[,-c(1:2)][,1:nPCs], 2,aggregate, by = list(KGped$SuperPopulation), median)
pop.medians<-cbind.data.frame(pop.medians)
rownames(pop.medians)<-pop.medians[,1]
pop.medians<-pop.medians[,seq(2,nPCs*2,2)]

## for each individual compare to each super population and find most similar
popDistsAll<-apply(pcas[,-c(1:2)][,1:nPCs], 1, calcPopDist, pop.medians)
popDistsAll<-t(popDistsAll)
predPop<-colnames(popDistsAll)[apply(popDistsAll, 1, which.min)]


ptCol<-rainbow(nSuperPops)[as.factor(predPop)]
pdf("PCAplotwith1KGpredictedPopulations.pdf", width = 10, height = 10)
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", type = "n")
points(pcas[which(ptType == 20),3], pcas[which(ptType == 20),4], pch = 20, col = ptCol[which(ptType == 20)])
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", type = "n")
points(pcas[which(ptType == 3),3], pcas[which(ptType == 3),4], pch = 3, col = ptCol[which(ptType == 3)])
dev.off()

## look at predications of our sample

outPred<-cbind(pcas[,1:2], predPop)
outPred<-outPred[which(ptType == 3),]
write.csv(table(outPred$predPop), "TablePredictedPopulations.csv")
write.csv(outPred, "PredictedPopulations.csv")

