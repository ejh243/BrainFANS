## calculate correlations between total and each cell fraction


source("rmdConfig.mrc")
setwd(dataDir)
load(normData)

## remove rs probes
celltypenormbeta<-celltypenormbeta[grep("rs", rownames(celltypenormbeta), invert = TRUE),]

## exclude IRF8+ve as not present in sufficient sample numbers
pheno<-pheno[pheno$Cell.type %in% c("Double -ve","NeuN +ve","Sox10 +ve", "Total"),]
celltypenormbeta<-celltypenormbeta[,pheno$Basename]
pheno$Cell.type<-factor(pheno$Cell.type)

## reformat as a list with matrix for each cell type, with samples in the same order
allSamples<-unique(pheno$Indidivual.ID)

betas.list<-list("Total" = celltypenormbeta[,which(pheno$Cell.type == "Total")[match(allSamples, pheno$Indidivual.ID[which(pheno$Cell.type == "Total")])]],
"NeuN" = celltypenormbeta[,which(pheno$Cell.type == "NeuN +ve")[match(allSamples, pheno$Indidivual.ID[which(pheno$Cell.type == "Total")])]],
"Sox10" = celltypenormbeta[,which(pheno$Cell.type == "Sox10 +ve")[match(allSamples, pheno$Indidivual.ID[which(pheno$Cell.type == "Total")])]],
"DNeg" = celltypenormbeta[,which(pheno$Cell.type == "Double -ve")[match(allSamples, pheno$Indidivual.ID[which(pheno$Cell.type == "Total")])]])

res<-matrix(data = NA, ncol = 3, nrow = nrow(celltypenormbeta))
colnames(res)<-c("NeuN", "Sox10", "DNeg")
rownames(res)<-rownames(celltypenormbeta)
for(i in 1:nrow(celltypenormbeta)){
	res[i,1]<-cor(betas.list[["Total"]][i,], betas.list[["NeuN"]][i,], use = "p")
	res[i,2]<-cor(betas.list[["Total"]][i,], betas.list[["Sox10"]][i,], use = "p")
	res[i,3]<-cor(betas.list[["Total"]][i,], betas.list[["DNeg"]][i,], use = "p")
}

save(res, file="Analysis/CellType/CTCorwithTotal.rdata")

