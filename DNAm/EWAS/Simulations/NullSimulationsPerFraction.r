## run null EWAS for each cell fraction to compare potential inflation effects
## limit data to same number of samples
## each each fraction only contains one individuals only need to use LM model

runEWAS<-function(row,pheno, status){

	modelLM<-lm(row ~ status + pheno$CCDNAmAge + pheno$Sex + pheno$Tissue.Centre)
	return(c(summary(modelLM)$coefficients["status1", 4]))
}


library(MASS)
#library(GenABEL)
library(doParallel)

nSim<-1000

args<-commandArgs(trailingOnly = TRUE)
source(args[1])

setwd(dataDir)

load(normData)

## remove cell types with less than 10 samples
nSample<-table(pheno$Cell.type)
pheno<-pheno[pheno$Cell.type %in% names(nSample[which(nSample > 9)]),]

celltypenormbeta<-celltypenormbeta[,pheno$Basename]

cellTypes<-unique(pheno$Cell.type)

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))


## first randomly assign cases control status, do we see inflation?
## fix so each cell fraction has same number of samples
minN<-min(table(as.character(pheno$Cell.type)))
allSim<-NULL
for(type in cellTypes){
	tmp.allSim<-NULL
	index<-which(pheno$Cell.type == type)
	index<-sample(index, minN)
	pheno.sub<-pheno[index,]
	celltypenormbeta.sub<-celltypenormbeta[,index]
	for(simNum in 1:nSim){
		status<-rep(0,nrow(pheno.sub))

		for(centre in unique(pheno.sub$Tissue.Centre)){
			ids<-unique(pheno.sub$Indidivual.ID[which(pheno.sub$Tissue.Centre == centre)])
			cases<-sample(ids, floor(length(ids)/2))
			status[pheno.sub$Indidivual.ID %in% cases]<-1
		}
		status<-as.factor(status)
		outtab<-foreach(i=1:nrow(celltypenormbeta.sub), .combine = "cbind") %dopar% runEWAS(celltypenormbeta.sub[i,], pheno.sub, status)
		outtab<-t(outtab)

		rownames(outtab)<-rownames(celltypenormbeta)
		tmp.allSim<-cbind(tmp.allSim, outtab)
	}
	
	allSim[[type]]<-tmp.allSim
}

save(allSim, file = "nullSimulationsPerFraction.rdata")
