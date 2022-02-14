## script to plot PCA output and colour by 1000G ethnicity

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


args<-commandArgs(TRUE)

setwd(args[1])
refFolder<-args[2]
subset<-args[3] ## use , to separate multiple files


if(!is.na(subset)){
	subset = unlist(strsplit(subset, ","))
	#print(subset)
	keepFiles<-sapply(subset, read.table, simplify=FALSE)
	
}


pcas<-read.table("mergedw1000G.pca.eigenvec", stringsAsFactors = FALSE)
KGped<-read.table(paste(refFolder, "/20130606_g1k.ped", sep = ""), stringsAsFactors = FALSE, header = TRUE, sep = "\t")
popInfo<-read.table(paste(refFolder, "/PopInfo.txt", sep = ""), stringsAsFactors = FALSE, header = TRUE, sep = "\t") ## table made from 1000G website

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

if(length(args) > 2){
	for(each in subset){
		par(mfrow = c(2,2))
		filterIndex<-which(pcas[,1] %in% keepFiles[[each]][,1] | !is.na(KGped[,1]))

		plot(pcas[filterIndex,3], pcas[filterIndex,4], xlab = "PC1", ylab = "PC2", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(pcas[filterIndex,3], pcas[filterIndex,5], xlab = "PC1", ylab = "PC3", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(pcas[filterIndex,3], pcas[filterIndex,6], xlab = "PC1", ylab = "PC4", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
		fileName<-tail(unlist(strsplit(each, "/")), n = 1)
		title(sub = paste(fileName, "\n(", nrow(keepFiles[[each]]), " individuals)", sep = ""))
		legend("center", pch = 16, col = rainbow(nPops), levels(as.factor(KGped$Population)), cex = 1.5, ncol = 3)
	}
}

## alternatively plot "super populations"


ptCol<-rainbow(nSuperPops)[as.factor(KGped$SuperPopulation)]
ptCol[is.na(ptCol)]<-"black"
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,5], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,3], pcas[,6], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)

if(length(args) > 2){
	for(each in subset){
		
		filterIndex<-which(pcas[,1] %in% keepFiles[[each]][,1] | !is.na(KGped[,1]))

		plot(pcas[filterIndex,3], pcas[filterIndex,4], xlab = "PC1", ylab = "PC2", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(pcas[filterIndex,3], pcas[filterIndex,5], xlab = "PC1", ylab = "PC3", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(pcas[filterIndex,3], pcas[filterIndex,6], xlab = "PC1", ylab = "PC4", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
		fileName<-tail(unlist(strsplit(each, "/")), n = 1)
		title(sub = paste(fileName, "\n(", nrow(keepFiles[[each]]), " individuals)", sep = ""))
		legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)
	}
}

dev.off()



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
plot(1:20,nMatches/sum(!is.na(KGped$SuperPopulation))*100, xlab = "nPCs", ylab = "Percentage Correct")
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

## calculate a quality score for prediction
## ideally want one population much closer than the others
rangeDist<-t(diff(apply(popDistsAll,1,range)))
qsPred<-(apply(popDistsAll,1,quantile, 0.25)-apply(popDistsAll,1,min))/rangeDist
pdf("BoxplotPrePopQCScores.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
boxplot(qsPred ~ KGped$SuperPopulation, col = rainbow(5), xlab = "Known populations")
boxplot(qsPred ~ predPop, col = rainbow(5), xlab = "Predicted populations")
dev.off()
## can define thresholds for each population based on 1000 genomes samples
#pop99<-aggregate(popDistsAll, by = list(KGped$SuperPopulation), quantile, 0.95)
#pop99Thres<-diag(as.matrix(pop99[,-1]))
#names(pop99Thres)<-colnames(popDistsAll)
#table(apply(popDistsAll, 1, min) < pop99Thres[predPop],as.factor(predPop))


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

if(length(args) > 2){
	for(each in subset){
		layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
		filterIndex<-which(pcas[,1] %in% keepFiles[[each]][,1] | !is.na(KGped[,1]))

		plot(pcas[filterIndex,3], pcas[filterIndex,4], xlab = "PC1", ylab = "PC2", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(pcas[filterIndex,3], pcas[filterIndex,5], xlab = "PC1", ylab = "PC3", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(pcas[filterIndex,3], pcas[filterIndex,6], xlab = "PC1", ylab = "PC4", pch = ptType[filterIndex], col = ptCol[filterIndex])
		plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
		fileName<-tail(unlist(strsplit(each, "/")), n = 1)
		title(sub = paste(fileName, "\n(", nrow(keepFiles[[each]]), " individuals)", sep = ""))
		#legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(KGped$SuperPopulation)), cex = 1.5)
		
		plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", type = "n")
		points(pcas[which(ptType == 20),3], pcas[which(ptType == 20),4], pch = 20, col = ptCol[which(ptType == 20)])
		plot(pcas[,3], pcas[,4], xlab = "PC1", ylab = "PC2", type = "n")
		points(pcas[which(pcas[,1] %in% keepFiles[[each]][,1]),3], pcas[which(pcas[,1] %in% keepFiles[[each]][,1]),4], pch = 3, col = ptCol[which(pcas[,1] %in% keepFiles[[each]][,1])])

	}
}

dev.off()


## look at predications of our sample
outPred<-cbind(pcas[,1:2], predPop, qsPred, popDistsAll)
outPred<-outPred[which(ptType == 3),]
write.csv(table(outPred$predPop), "TablePredictedPopulations.csv")
write.csv(outPred, "PredictedPopulations.csv", quote = FALSE, row.names = FALSE)

if(length(args) > 2){
	for(each in subset){
		fileName<-tail(unlist(strsplit(each, "/")), n = 1)
		outPred.tmp<-outPred[which(outPred[,1] %in% keepFiles[[each]][,1]),]
		write.csv(table(outPred.tmp$predPop), paste("TablePredictedPopulations", fileName, ".csv", sep = ""))
	}
}
		
