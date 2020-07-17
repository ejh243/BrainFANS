## for each cell type (CT) calculate variance characteristic score
## first adjust betas for differences in mean level of DNAm
## second use levene's test to compare variance in one cell type against all others
## one sided test 


library(car)
library(parallel)
cl <- makeCluster(16)
clusterExport(cl,"leveneTest")

adjustCT<-function(row, sampleType){
	return(residuals(lm(row ~ sampleType)))
}

varCharScore<-function(row, indVar){
	return(leveneTest(row ~ indVar, center = "median")[1,3])
}


source("rmdConfig.mrc")
setwd(dataDir)
load(normData)

## remove rs probes
celltypenormbeta<-celltypenormbeta[grep("rs", rownames(celltypenormbeta), invert = TRUE),]

## exclude IRF8+ve as not present in sufficient sample numbers
pheno<-pheno[pheno$Cell.type %in% c("Double -ve","NeuN +ve","Sox10 +ve"),]
celltypenormbeta<-celltypenormbeta[,pheno$Basename]
pheno$Cell.type<-factor(pheno$Cell.type)


betasAdj<-t(apply(celltypenormbeta, 1, adjustCT, pheno$Cell.type))

pheno$Cell.type<-as.factor(as.character(pheno$Cell.type))

out<-matrix(data = NA, nrow = nrow(betasAdj), ncol = 3)
colnames(out)<-levels(pheno$Cell.type)
rownames(out)<-rownames(celltypenormbeta)
for(each in levels(pheno$Cell.type)){
	indVar<-rep(0, ncol(betasAdj))
	indVar[which(pheno$Cell.type == each)]<-1
	indVar<-as.factor(indVar)
	tmp<-parApply(cl, betasAdj, 1, varCharScore, indVar)
	out[,each]<-tmp
	out[,each]<--log10(out[,each])
	## as non-parametric test can test for spcific direction so instead will use sign to indicate direction of effect
	## calculate difference in sd to determine sign
	coef<-sign(parApply(cl, betasAdj[,which(indVar == 1)], 1, sd)-parApply(cl,betasAdj[,which(indVar == 0)], 1, sd))
	out[,each]<-out[,each]*coef
}

save(out, file = "Analysis/CellType/VarianceCharacteristicScores.rdata")


