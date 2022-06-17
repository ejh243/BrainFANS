## run this script to confirm sorted cell fractions are labelled correctly and data clusters as expected
## use QC metrics for phenotype data
## within adult samples do PCA to cluster sample types


args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
config <- args[2]

source(config)

gdsFile <-paste0(dataDir, "/2_gds/raw.gds")
qcOutFolder<-paste0(dataDir, "/2_gds/QCmetrics")
qcData <-paste0(dataDir, "/2_gds/QCmetrics/QCmetrics.rdata")
genoFile <- paste0(dataDir, "/0_metadata/epicSNPs.raw")


library(bigmelon)


setwd(dataDir)


gfile<-openfn.gds(gdsFile, readonly = FALSE)

## filter samples
QCSum<-read.csv(paste0(dataDir, "/2_gds/QCmetrics/PassQCStatusAllSamples.csv"), row.names = 1, stringsAsFactors = FALSE)
passQC<-QCSum$Basename[QCSum[,"passQCS2"]]

load(qcData)
QCmetrics.all<-QCmetrics
QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

rawbetas<-gfile[,, node = "betas"]
rawbetas<-rawbetas[,match(passQC, colnames(rawbetas))]
auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
rawbetas<-rawbetas[auto.probes,]

cellTypes<-unique(QCmetrics$Cell.type)
cellTypes<-cellTypes[!is.na(cellTypes)]
cellTypes<-sort(cellTypes)



## filter out NAs
rawbetas<-na.omit(rawbetas)

## filter out SNPs
rawbetas<-rawbetas[-grep("rs", rownames(rawbetas)),]

print("Calculating PCs from beta matrix")

## Run PCA across all cell types
pca <- prcomp(t(rawbetas))
betas.scores = pca$x
colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
betas.pca<-pca$sdev^2/sum(pca$sdev^2)
betas.scores<-betas.scores[,which(betas.pca > 0.01)]

## Following methods rely on the ability to define an average profile of a cell type so need to check if any outliers first
print("Identifying outliers")

## calculate Studentized scores for each sample compared to it's cell type mean
studentPCA<-list()
studentPCA[[1]]<-matrix(NA, nrow = nrow(QCmetrics), ncol = ncol(betas.scores))
for(i in 1:nrow(QCmetrics)){
	keep<-QCmetrics$Basename[which(QCmetrics$Cell.type == QCmetrics$Cell.type[i])]
	## exclude itself
	keep<-keep[keep != QCmetrics$Basename[i]]
	cell.means<-colMeans(betas.scores[keep,])
	cell.sd<-apply(betas.scores[keep,], 2, sd)
	studentPCA[[1]][i,]<-c((betas.scores[i,]-cell.means)/cell.sd)
}


## identify outliers from first 2 PCs only
studentCTOutlier<-rowMax(abs(studentPCA[[1]][,1:2])) > studentThres
QCmetrics$maxStudentPCA<-rowMax(abs(studentPCA[[1]][,1:2]))
QCmetrics$studentCTOutlier<-studentCTOutlier


## for each sample quantify how similar to cell type average it is

cellMeanPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier,], by = list(QCmetrics$Cell.type[!QCmetrics$studentCTOutlier]), mean, na.rm = TRUE)
rownames(cellMeanPCA)<-cellMeanPCA[,1]
cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))

cellSDPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier,], by = list(QCmetrics$Cell.type[!QCmetrics$studentCTOutlier]), sd, na.rm = TRUE)
rownames(cellSDPCA)<-cellSDPCA[,1]
cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))

maxSD<-rep(NA, nrow(QCmetrics))
for(i in 1:nrow(QCmetrics)){
  maxSD[i]<-max(abs(betas.scores[i,1:2]-cellMeanPCA[1:2,QCmetrics$Cell.type[i]])/cellSDPCA[1:2,QCmetrics$Cell.type[i]])

}


## For each individual calculate a FACs efficency score
uniqueIDs<-unique(QCmetrics[,c("Indidivual.ID", "Sex", "Age", "Phenotype", "Tissue.Centre")])
indFACSEff<-indFACSEff<-aggregate(maxSD[which(QCmetrics$Cell.type != "Total")], by = list(QCmetrics$Indidivual.ID[which(QCmetrics$Cell.type != "Total")]), FUN = median)
nFACs<-table(QCmetrics$Indidivual.ID[QCmetrics$Cell.type != "Total"])

uniqueIDs<-cbind(uniqueIDs, indFACSEff$x[match(uniqueIDs$Indidivual.ID, as.character(indFACSEff$Group.1))], as.numeric(nFACs[uniqueIDs$Indidivual.ID]))
colnames(uniqueIDs)<-c("Indidivual.ID", "Sex", "Age", "Phenotype", "Tissue.Centre", "FACsEffiency", "nFACS")

write.csv(uniqueIDs, paste0(qcOutFolder, "/IndividualFACsEffciencyScores.csv"))

## exclude individuals who very poor FACS sorts

QCmetrics$passFACS<-QCmetrics$Indidivual.ID %in% uniqueIDs$Indidivual.ID[which(uniqueIDs$FACsEffiency < 5)]

## repeat outlier identification
## calculate Studentized scores for each sample compared to it's cell type mean
studentPCA[[2]]<-matrix(NA, nrow = nrow(QCmetrics), ncol = ncol(betas.scores))
for(i in 1:nrow(QCmetrics)){
	if(QCmetrics$passFACS[i]){
		keep<-QCmetrics$Basename[which(QCmetrics$Cell.type == QCmetrics$Cell.type[i] & QCmetrics$passFACS)]
		## exclude itself
		keep<-keep[keep != QCmetrics$Basename[i]]
		cell.means<-colMeans(betas.scores[keep,])
		cell.sd<-apply(betas.scores[keep,], 2, sd)
		studentPCA[[2]][i,]<-c((betas.scores[i,]-cell.means)/cell.sd)
	}
}

## identify outliers from first 2 PCs only
studentPCA[[2]][is.na(studentPCA[[2]])]<-100
studentCTOutlier<-rowMax(abs(studentPCA[[2]][,1:2])) > studentThres
QCmetrics$maxStudentPCA2<-rowMax(abs(studentPCA[[2]][,1:2]))
QCmetrics$maxStudentPCA2[which(QCmetrics$maxStudentPCA2 == 100)]<-NA
QCmetrics$studentCTOutlier2<-studentCTOutlier


## for each sample quantify how similar to the cell type average it is
cellMeanPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier2,], by = list(QCmetrics$Cell.type[!QCmetrics$studentCTOutlier2]), mean, na.rm = TRUE)
rownames(cellMeanPCA)<-cellMeanPCA[,1]
cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))

cellSDPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier2,], by = list(QCmetrics$Cell.type[!QCmetrics$studentCTOutlier2]), sd, na.rm = TRUE)
rownames(cellSDPCA)<-cellSDPCA[,1]
cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))

maxSD<-rep(NA, nrow(QCmetrics))
for(i in 1:nrow(QCmetrics)){
  maxSD[i]<-max(abs(betas.scores[i,1:2]-cellMeanPCA[1:2,QCmetrics$Cell.type[i]])/cellSDPCA[1:2,QCmetrics$Cell.type[i]])

}

QCmetrics<-cbind(QCmetrics, maxSD)

## use Mahalanobis to calculate distance with cell type medians from PCAs 
## exclude outliers identified above from calculation of median profile

print("Calculating Mahalanobis distance")

cellMedPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier2,], by = list(QCmetrics$Cell.type[!QCmetrics$studentCTOutlier2]), median, na.rm = TRUE)
rownames(cellMedPCA)<-cellMedPCA[,1]
cellMedPCA<-as.matrix(t(cellMedPCA[,-1]))

## calc cell type covariance matrices
cov_sigmaPCA<-list()
for(each in colnames(cellMedPCA)){
	## only possible if > 1 sample
	if(sum(QCmetrics$Cell.type == each) > 1){
		cov_sigmaPCA[[each]]<-cov(betas.scores[which(QCmetrics$Cell.type == each & !QCmetrics$studentCTOutlier2),1:2], use = "p")
	}
}

mahDistPCA<-matrix(data = 10^9, ncol = ncol(cellMedPCA), nrow = nrow(QCmetrics))
colnames(mahDistPCA)<-colnames(cellMedPCA)
for(each in colnames(cellMedPCA)){
	if(!is.null(cov_sigmaPCA[[each]]) && det(cov_sigmaPCA[[each]]) > 1){
		mahDistPCA[,each]<-mahalanobis(betas.scores[,1:2], cellMedPCA[,each], cov_sigmaPCA[[each]], na.rm = TRUE,tol=1e-20)
	}
}

closestCellTypePCA<-colnames(mahDistPCA)[unlist(apply(mahDistPCA, 1, which.min))]
closestLabelledCellType<-QCmetrics$Cell.type == closestCellTypePCA

## As two different antibodies used for neuronal cells used, accept if either predicted as the other
neunIndex<-which(QCmetrics$Cell.type %in% neunCT)
closestLabelledCellType[neunIndex[closestCellTypePCA[neunIndex] %in% neunCT]]<-TRUE

## SATB2- are a composition of non neuronal cells; accept so long as predicted as non-neuronal
satb2NegIndex<-which(QCmetrics$Cell.type == "SATB2-")
closestLabelledCellType[satb2NegIndex[!closestCellTypePCA[satb2NegIndex] %in% neunCT]]<-TRUE

QCmetrics<-cbind(QCmetrics, closestCellTypePCA, closestLabelledCellType)
write.csv(QCmetrics[which(closestLabelledCellType == "FALSE"),], paste0(qcOutFolder, "/SamplesPredictedDiffCellTypePCAMahDist.csv"))


## additionally check sample overlaps with mean profile of that cell type
## again exclude outliers identified above from calculating the average profile
print("Confirming sample overlap with group polytope")
outlierThres<-c(1.5,2,2.5,3)
pcaClassify<-list("predictCellType" = matrix(data = NA, nrow = nrow(QCmetrics), ncol = length(outlierThres)), "withinSDMean" = matrix(data = NA, nrow = nrow(QCmetrics), ncol = length(outlierThres)))
for(thres in outlierThres){
	studentCTOutlier<-QCmetrics$maxStudentPCA2 > thres

	cellMeanPCA<-aggregate(betas.scores[!studentCTOutlier,], by = list(QCmetrics$Cell.type[!studentCTOutlier]), mean, na.rm = TRUE)
	rownames(cellMeanPCA)<-cellMeanPCA[,1]
	cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))

	cellSDPCA<-aggregate(betas.scores[!studentCTOutlier,], by = list(QCmetrics$Cell.type[!studentCTOutlier]), sd, na.rm = TRUE)
	rownames(cellSDPCA)<-cellSDPCA[,1]
	cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))

	lowerBound<-cellMeanPCA-nSDThres*cellSDPCA
	upperBound<-cellMeanPCA+nSDThres*cellSDPCA
	

	# assumes cell types are distinct so just consider set of non-overlappping cell types
	pcaClassifyDistinctCT<-rep("", nrow(QCmetrics))
	for(each in predDistinctCT){
		indexClassify<-which(betas.scores[,1] < upperBound[1,each] & betas.scores[,1] > lowerBound[1,each] & betas.scores[,2] < upperBound[2,each] & betas.scores[,2] > lowerBound[2,each])
		pcaClassifyDistinctCT[indexClassify]<-paste(pcaClassifyDistinctCT[indexClassify], each, sep = ";")
	}

	pcaClassify[["predictCellType"]][,match(thres, outlierThres)]<-pcaClassifyDistinctCT

	## can overlap with multiple CT so check if it if within 2SD of mean of it's labelled CT
	withinSDMean<-rep(FALSE, nrow(QCmetrics))
	for(i in 1:nrow(QCmetrics)){
		if(QCmetrics$Cell.type[i] %in% colnames(upperBound)){
			withinSDMean[i]<-betas.scores[i,1] < upperBound[1,QCmetrics$Cell.type[i]] & betas.scores[i,1] > lowerBound[1,QCmetrics$Cell.type[i]] & betas.scores[i,2] < upperBound[2,QCmetrics$Cell.type[i]] & betas.scores[i,2] > lowerBound[2,QCmetrics$Cell.type[i]]
		}
	}
	
	pcaClassify[["withinSDMean"]][,match(thres, outlierThres)]<-withinSDMean
}

## pull out results for specified threshold
withinSDMean<-pcaClassify[["withinSDMean"]][,match(studentThres, outlierThres)]
pcaClassifyDistinctCT<-pcaClassify[["predictCellType"]][,match(studentThres, outlierThres)]

write.csv(QCmetrics[which(withinSDMean == "FALSE"),],paste0(qcOutFolder, "/SamplesPCAOutlierFromCellType.csv"))

## exclude individuals with sub optimal FACs
pcaClassifyDistinctCT[!QCmetrics$passFACS]<-FALSE
withinSDMean[!QCmetrics$passFACS]<-FALSE
QCmetrics<-cbind(QCmetrics, pcaClassifyDistinctCT,withinSDMean)

## keep all TOTAL samples
closestLabelledCellType[which(QCmetrics$Cell.type == "Total")]<-TRUE
QCmetrics$predLabelledCellType<-closestLabelledCellType
withinSDMean[which(QCmetrics$Cell.type == "Total")]<-TRUE
QCmetrics$withinSDMean<-withinSDMean

passCTCheck<-withinSDMean
passCTCheck[which(QCmetrics$Cell.type == "Total")]<-TRUE
QCmetrics$passCTCheck<-passCTCheck

## for time being take more stringent approach.
#add.gdsn(gfile, 'QCdata', val = QCmetrics, replace = TRUE)

## need to close gds file in order to open in another R session
closefn.gds(gfile)

save(betas.scores, mahDistPCA, pcaClassify, studentPCA, file = paste0(qcOutFolder,"/PCAAcrossAllCellTypes.rdata"))

print("Updating QC summary")

## add outcome to qc summary
QCmetrics<-QCmetrics[match(QCSum$Basename, QCmetrics$Basename),]

## add in data for samples excluded thus far
QCmetrics[which(is.na(QCmetrics$withinSDMean)),colnames(QCmetrics.all)]<-QCmetrics.all[which(is.na(QCmetrics$withinSDMean)),]

#QCSum<-QCSum[,-ncol(QCSum)]
QCSum<-cbind(QCSum, QCmetrics$passFACS, QCmetrics$passCTCheck,  QCmetrics$passFACS & QCmetrics$passCTCheck)
colnames(QCSum)[(ncol(QCSum)-2):ncol(QCSum)]<-c("passCTCheck", "passFACS", "passQCS3")

write.csv(QCSum, paste0(qcOutFolder,"/passQCStatusStage3AllSamples.csv"))
write.csv(QCmetrics, paste0(qcOutFolder,"/QCMetricsPostCellTypeClustering.csv"))