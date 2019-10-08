## quick QC to check array quality


setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/")
library(bigmelon)
library(gplots)

## load data

sampleSheet<-read.csv("SampleSheets/Run1June2019.csv")
sampleSheet$Basename<-paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
setwd("iDats")

for(each in sampleSheet$Basename){
	gfile <- iadd2(".", gds = 'run1.gds') ## this only reads in a single gds
}

 gfile <- openfn.gds('run1.gds')

### three quick checks:

msetEPIC <- readEPIC(idatPath="iDats", barcodes=sampleSheet$Basename, parallel = FALSE, force=T)

m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)
QCmetrics<-cbind(pheno,M.median, U.median)

plot(M.median, U.median, col = as.factor(sampleSheet$Cell.type), pch = 16)


betas.tmp<-betas(msetEPIC)
sigma<-apply(betas.tmp, 1, sd, na.rm = TRUE)
heatmap.2(betas.tmp[which(sigma > quantile(sigma, 0.995, na.rm = TRUE)),], trace = "none", dendrogram = "column", ColSideColors = as.factor(sampleSheet$Cell.type), labRow = "", labCol = "")

grep("rs", rownames(betas.tmp))
heatmap.2(betas.tmp[grep("rs", rownames(betas.tmp)),], trace = "none", dendrogram = "column", ColSideColors = c("black", "blue", "red", "forestgreen")[as.factor(sampleSheet$Cell.type)], labRow = "", labCol = sampleSheet$Individual)

Individual


## see if samples cluster by type
betas.tmp<-read.gdsn(betas(gfile))
betas.tmp<-read.gdsn(fData(gfile))
sigma<-apply(betas.tmp, 1, sd, na.rm = TRUE)

heatmap.2(betas.tmp[which(sigma > quantile(sigma, 0.95, na.rm = TRUE)),], trace = "none")
hclust(dist(betas.tmp[which(sigma > quantile(sigma, 0.99, na.rm = TRUE)),]))

## see if samples cluster by individual across snp probes

grep("rs", fData(gfile))