## calculate variance explained in total by all three cell fractions

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

## take total as outcome then all other cell types as variables

out<-matrix(data = NA, ncol = 2, nrow = nrow(betas.list[["Total"]]))
rownames(out)<-rownames(betas.list[["Total"]])
colnames(out)<-c("MultipleR2", "AdjustedR2")

for(i in 1:nrow(out)){
	model<-lm(betas.list[["Total"]][i,] ~ betas.list[["NeuN"]][i,] + betas.list[["Sox10"]][i,] + betas.list[["DNeg"]][i,])
	out[i,2]<-summary(model)$adj.r.squared
	out[i,1]<-summary(model)$r.squared
}

save(out, file="Analysis/CellType/TotalVarExplainAllCT.rdata")