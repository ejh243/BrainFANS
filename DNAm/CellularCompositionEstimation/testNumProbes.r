## script to develop neural cell type deconvolution algorithm
## uses houseman method with new ref data for purified brain cell populations
## editted original code from minfi to take a matrix rather than RGset - this means data is not preprocessed together.
## this script explores how number of probes influences the result
## uses two testing scenarios:
## 1. ability to predict purified cell type, 
## 2. ability to accurately estimate cellular composition of computationally simulated bulk profiles. 


classifyCellType<-function(row,thres, names){
	if(max(row) > thres){
		return(names[which(row > thres)])
	} else{
		return(NA)
	}
}

source("FunctionsForBrainCellProportionsPrediction.r")
library(minfi)
library(genefilter)
library(rafalib)
library(RColorBrewer)
source("rmdConfig.bdr")
setwd(dataDir) 
load(normData)


pheno$Cell.type<-gsub(" \\+ve", "", pheno$Cell.type) ## need to remove the "+" and "-"
pheno$Cell.type<-gsub(" -ve", "", pheno$Cell.type) ## need to remove the "+" and "-"

## parameters
probeSelect = "any" ## options "both" (select equal number of probes associated with hyper and hypo methylation) or "any" (just select probes based on significance regardless of direction)
cellTypes = c("Double","NeuN", "Sox10", "IRF8")
cellTypes<-sort(cellTypes)


## filter data to selected cellTypes
pheno<-pheno[which(pheno$Cell.type %in% cellTypes),]
celltypenormbeta<-celltypenormbeta[,pheno$Basename]

cellInd<-as.factor(pheno$Cell.type)

### Test 1: split into train and test 75% vs 25%
## compare performance across range of different number of probes
## outcome: correct prediction (classed at > 95% proportion).

nSamples<-round(table(pheno$Cell.type)*0.75)

rangeProbes<-c(seq(5,100, 5),seq(120,1000,20))
nSims<-10
nCorrect<-lapply(1:5, matrix,data = NA, ncol = nSims, nrow = length(rangeProbes))
names(nCorrect)<-c("All", cellTypes)
meanPredCorrect<-lapply(1:5, matrix,data = NA, ncol = nSims, nrow = length(rangeProbes))
names(meanPredCorrect)<-c("All", cellTypes)
nProbes<-matrix(data = NA, ncol = nSims, nrow = length(rangeProbes))

for(j in 1: nSims){
samples.train<-NULL
for(i in 1:length(nSamples)){
	samples.train<-c(samples.train, sample(which(pheno$Cell.type == names(nSamples)[i]),as.numeric(nSamples[i])))
}	
samples.test<-c(1:nrow(pheno))[!1:nrow(pheno) %in% samples.train]

betas.train<-celltypenormbeta[,pheno$Basename[samples.train]]
betas.test<-celltypenormbeta[,pheno$Basename[samples.test]]

cellInd.train<-cellInd[samples.train]
cellInd.test<-cellInd[samples.test]
		
## perform F-test to identify probes with differences between cell types.
## NB only need to do F-test once for each simulation
ffComp <- rowFtests(betas.train, cellInd.train)
prof <- vapply(
        X = splitit(cellInd.train),
        FUN = function(j) rowMeans2(betas.train, cols = j),
        FUN.VALUE = numeric(nrow(betas.train)))
r <- rowRanges(betas.train)
compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
names(compTable)[1] <- "Fstat"
names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
        c("low", "high", "range")
tIndexes <- splitit(cellInd.train)
tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(betas.train))
        x[i] <- 1
        return(rowttests(betas.train, factor(x)))
})

form <- as.formula(
       sprintf("y ~ %s - 1", paste(levels(cellInd.train), collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~ cellInd.train - 1))
colnames(phenoDF) <- sub("cellInd.train", "", colnames(phenoDF))

## select probes for deconvolution
## look at effect of number of probes
rowNum<-1
for(numProbes in rangeProbes){

	if (probeSelect == "any") {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-4, ]
            yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
            c(rownames(yAny)[seq(numProbes * 2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-4, ]
            yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
            yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
            c(rownames(yUp)[seq_len(numProbes)],
              rownames(yDown)[seq_len(numProbes)])
        })
    }

	## if insufficient number of probes to select from NAs creep into probeList
	probeList<-lapply(probeList, na.omit)
    trainingProbes <- unique(unlist(probeList))
    p <- betas.train[trainingProbes,]
	
	# > 2 groups solution
    tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts

	nProbes[rowNum,j]<-length(trainingProbes)
	
	## generate prediction in test data
	counts <- projectCellType(betas.test[rownames(coefEsts), ], coefEsts)
	predCT<-apply(counts, 1, classifyCellType, 0.9, colnames(counts))
	nCorrect[["All"]][rowNum,j]<-sum(predCT == cellInd.test, na.rm = TRUE)/length(samples.test)
	propCorrectCT<-diag(apply(counts, 1, "[", match(cellInd.test, colnames(counts))))
	meanPredCorrect[["All"]][rowNum,j]<-mean(propCorrectCT, na.rm = TRUE)
	for(each in cellTypes){
		index<-which(cellInd.test == each)
		nCorrect[[each]][rowNum,j]<-sum((predCT == cellInd.test)[index], na.rm = TRUE)/length(index)
		meanPredCorrect[[each]][rowNum,j]<-mean(propCorrectCT[index], na.rm = TRUE)
	}
	rowNum<-rowNum+1
}

}
plotCols<-brewer.pal(length(cellTypes)+1, "Set1")
names(plotCols)<-c("All", cellTypes)
## calculate median and quartiles across simulations

sumStats.1<-lapply(nCorrect, apply, 1, quantile, c(0.25,0.5,0.75))
sumStats.2<-lapply(meanPredCorrect, apply, 1, quantile, c(0.25,0.5,0.75))
pdf("TestCase1SummaryPlots.pdf", width = 8, height = 8)
par(mar = c(4,4,1,1))
par(mfrow = c(2,2))
y_lim<-range(unlist(lapply(sumStats.1, range)))
plot(nProbes[,1], t(sumStats.1[["All"]])[,2], type = "l", ylim = y_lim, xlab = "nProbes", ylab = "Proportion Correct", col = plotCols[1], lwd = 2)
for(each in cellTypes){
	lines(nProbes[,1], t(sumStats.1[[each]])[,2], col = plotCols[each], lwd = 2)
}
#legend("bottomleft", names(plotCols), lty = 1, col = plotCols)

y_lim<-range(unlist(lapply(sumStats.2, range)))
y_lim[1]<-y_lim[1]-diff(y_lim)*0.2
plot(nProbes[,1], t(sumStats.2[["All"]])[,2], type = "l", ylim = y_lim, xlab = "nProbes", ylab = "Mean predicted proportion", col = plotCols[1], lwd = 2)
for(each in cellTypes){
	lines(nProbes[,1], t(sumStats.2[[each]])[,2], col = plotCols[each], lwd = 2)
}
legend("bottomleft", names(plotCols), lty = 1, col = plotCols)

y_lim<-range(unlist(lapply(sumStats.1, range)))
plot(nProbes[,1], t(sumStats.1[["All"]])[,2], type = "l", ylim = y_lim, xlab = "nProbes", ylab = "Proportion Correct", col = plotCols[1], lwd = 2, xlim = c(0,2500))
for(each in cellTypes){
	lines(nProbes[,1], t(sumStats.1[[each]])[,2], col = plotCols[each], lwd = 2)
}
#legend("bottomleft", names(plotCols), lty = 1, col = plotCols)

y_lim<-range(unlist(lapply(sumStats.2, range)))
plot(nProbes[,1], t(sumStats.2[["All"]])[,2], type = "l", ylim = y_lim, xlab = "nProbes", ylab = "Mean predicted proportion", col = plotCols[1], lwd = 2, xlim = c(0,2500))
for(each in cellTypes){
	lines(nProbes[,1], t(sumStats.2[[each]])[,2], col = plotCols[each], lwd = 2)
}
dev.off()

### Test 2: select individuals with all 4 fractions, to compute bulk simulations across full full matrix of possible proportions
## train in remaining samples

nSamp<-table(pheno$Indidivual.ID)
compIDs<-names(nSamp)[which(nSamp == length(cellTypes))]

# create a matrix fo all combinations for three cell types
# keep those that total 1 or less
# add a 4th cell type that makes up the difference to one. 
intervals<-seq(0,1,0.05)
propProfile<-expand.grid(rep(list(intervals), 3))
tot<-rowSums(propProfile)
propProfile<-propProfile[which(tot < 1),]
leftOver<-1-rowSums(propProfile)
propProfile<-cbind(propProfile, leftOver)
colnames(propProfile)<-cellTypes

bulkSims<-data.frame(matrix(vector(), 0, 2+(2*length(cellTypes)), dimnames=list(c(), c("ID","NumProbes", paste0("ActualProp", cellTypes),paste0("EstimatedProp", cellTypes)))),stringsAsFactors=F)

for(ID in compIDs){

	pureProfiles<-celltypenormbeta[, which(pheno$Indidivual.ID == ID)]
	colnames(pureProfiles)<-cellInd[which(pheno$Indidivual.ID == ID)]
	pureProfiles<-pureProfiles[,cellTypes]
	
	## create a matrix of bulk profile weighted by proportions define earlier
	bulkProfiles <- pureProfiles %*% t(propProfile)
	
	## train model in all other samples
	betas.train<-celltypenormbeta[,which(pheno$Indidivual.ID != ID)]
	cellInd.train<-cellInd[which(pheno$Indidivual.ID != ID)]
	## perform F-test to identify probes with differences between cell types.
	## NB only need to do F-test once for each simulation
	ffComp <- rowFtests(betas.train, cellInd.train)
	prof <- vapply(
			X = splitit(cellInd.train),
			FUN = function(j) rowMeans2(betas.train, cols = j),
			FUN.VALUE = numeric(nrow(betas.train)))
	r <- rowRanges(betas.train)
	compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
	names(compTable)[1] <- "Fstat"
	names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
			c("low", "high", "range")
	tIndexes <- splitit(cellInd.train)
	tstatList <- lapply(tIndexes, function(i) {
			x <- rep(0,ncol(betas.train))
			x[i] <- 1
			return(rowttests(betas.train, factor(x)))
	})

	form <- as.formula(
		   sprintf("y ~ %s - 1", paste(levels(cellInd.train), collapse = "+")))
	phenoDF <- as.data.frame(model.matrix(~ cellInd.train - 1))
	colnames(phenoDF) <- sub("cellInd.train", "", colnames(phenoDF))

	## select probes for deconvolution
	## look at effect of number of probes
	rowNum<-1
	for(numProbes in rangeProbes){

		if (probeSelect == "any") {
			probeList <- lapply(tstatList, function(x) {
				y <- x[x[, "p.value"] < 1e-4, ]
				yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
				c(rownames(yAny)[seq(numProbes * 2)])
			})
		} else {
			probeList <- lapply(tstatList, function(x) {
				y <- x[x[, "p.value"] < 1e-4, ]
				yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
				yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
				c(rownames(yUp)[seq_len(numProbes)],
				  rownames(yDown)[seq_len(numProbes)])
			})
		}

		## if insufficient number of probes to select from NAs creep into probeList
		probeList<-lapply(probeList, na.omit)
		trainingProbes <- unique(unlist(probeList))
		p <- betas.train[trainingProbes,]
		
		# > 2 groups solution
		tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
		coefEsts <- tmp$coefEsts

		nProbes[rowNum,j]<-length(trainingProbes)
		
		## generate prediction in test data
		counts <- projectCellType(bulkProfiles[rownames(coefEsts), ], coefEsts)
		tmp<-cbind(ID, length(trainingProbes), propProfile, counts)
		colnames(tmp)<-c("ID","NumProbes", paste0("ActualProp", cellTypes),paste0("EstimatedProp", cellTypes))
		bulkSims<-rbind(bulkSims,tmp)
	}
}

## need to fix precision error
bulkSims[,grep("ActualProp", colnames(bulkSims))]<-signif(bulkSims[,grep("ActualProp", colnames(bulkSims))],2)

## calculate RMSE
diffTab<-abs(bulkSims[,grep("ActualProp", colnames(bulkSims))] - bulkSims[,grep("EstimatedProp", colnames(bulkSims))])
colnames(diffTab)<-paste0("Diff", cellTypes)
bulkSims<-cbind(bulkSims, diffTab)
plotValues<-aggregate(as.formula(sprintf("cbind(%s) ~ NumProbes", paste(paste0("Diff", cellTypes), collapse = ","))), bulkSims,FUN = "median")

## focus on numProbes ~ 2000
## plot as a function of actual proportion

bulkSims.sub<-bulkSims[which(bulkSims$NumProbes <= 2000 & bulkSims$NumProbes >= 1900),]
medErrorByProp<-list()
for(each in cellTypes){
	medErrorByProp[[each]]<-aggregate(as.formula(sprintf("%s ~ %s", paste0("Diff", each), paste0("ActualProp", each))), bulkSims.sub,FUN = "median")
}


pdf("TestCase2SummaryPlots.pdf", width = 8, height = 8)
par(mar = c(4,4,1,1))
par(mfrow = c(2,1))
y_lim<-range(unlist(lapply(plotValues[,-1], range)))
plot(plotValues[,1], plotValues[,2], type = "n", ylim = y_lim, xlab = "nProbes", ylab = "Median Error")
for(i in 2:5){
	lines(plotValues[,1], plotValues[,i], col = plotCols[i], lwd = 2)
}

y_lim<-c(0, max(unlist(lapply(lapply(medErrorByProp, "[",2), max))))
plot(medErrorByProp[[1]][,1], medErrorByProp[[1]][,2], type = "n", ylim = y_lim, xlab = "Actual Proportion", ylab = "Median Error")
for(each in cellTypes){
	lines(medErrorByProp[[each]][,1], medErrorByProp[[each]][,2], col = plotCols[each], lwd = 2)
}
y_lim[2]<-y_lim[2]+diff(y_lim)*0.2

legend("topleft", names(plotCols), lty = 1, col = plotCols)

dev.off()

save(bulkSims, nCorrect, nProbes, meanPredCorrect, file = "SimulationsNumProbes.rdata")