##---------------------------------------------------------------------#
##
## Title: Look for overlap between EWAs and GWAS results
##
## Purpose of script: 
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refPath <- args[2]
gwasFile <- args[3]

normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resPath <- file.path(dataDir, "4_analysis/EWAS")

cellTypes <- c("Double-", "NeuN+", "Sox10+")

#----------------------------------------------------------------------#
# DEFINE FUNCTIONS
#----------------------------------------------------------------------#

calcCovMatrix<-function(x){
   ### x is a vector of correlations
   y<-vector(length = length(x))
   for(i in 1:length(x)){
      if(x[i] <= 1 & x[i] >= 0){
         y[i]<-x[i]*(3.25 + 0.75*x[i])
    } else {
      if(x[i] <= 0 & x[i] >= -0.5){
         y[i]<-x[i]*(3.27 + 0.71*x[i])
      }
    }
   
  }
  return(y)
}

brownsP<-function(covar, pval){
	## covar is vector of covariances between all pairs of tests NOTE not correlations
	## pval is vector of p values
	ntests<-length(pval)
	var_xsq<-sum(covar)*2 + 4*ntests 	# equation (3)
	exp_xsq<-2*ntests	#equation (2)

	## estimate parameters for chi square distribution using Brown's notation
	f = (2*(exp_xsq^2))/var_xsq
	c = var_xsq/(2*exp_xsq)

	##### NOTE: doesn't match Brown's answer but matches my answer by hand
	chi_sq<- -2*sum(log(pval))

	### to obtain p value
	test.stat<-chi_sq/c
	browns.pval<-1-pchisq(test.stat, df = f)
	return(c(browns.pval, test.stat, f))
}

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(GenomicRanges)
library(ggplot2)
library(ggpubr)

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#


res<-list()
for (each in cellTypes) {
    load(file.path(resPath, paste0(each, "LM.rdata")))
    res[[each]] <- outtab
}
rm(outtab)


setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

sampleNumbers<-table(QCmetrics$Phenotype, QCmetrics$Cell.type)


#----------------------------------------------------------------------#
# REMOVE CROSS HYB & SNP PROBES
#----------------------------------------------------------------------#


crosshyb <- read.table(file.path(refPath, "CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)
tofilter <- read.csv(file.path(refPath, "EPICArrayProbesToFilter.csv"), stringsAsFactors = FALSE)
snpProbes <- read.table(file.path(refPath, "SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
crosshyb2 <- read.csv(file.path(refPath, "Pidsley_SM1.csv"), stringsAsFactors = FALSE)
snpProbes2 <- read.csv(file.path(refPath, "Pidsley_SM4.csv"), stringsAsFactors = FALSE)
snpProbes3 <- read.csv(file.path(refPath, "Pidsley_SM5.csv"), stringsAsFactors = FALSE)
snpProbes4 <- read.csv(file.path(refPath, "Pidsley_SM6.csv"), stringsAsFactors = FALSE)

snpProbes <- snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$AF > 0.01), ]
snpProbes2 <- snpProbes2[which(snpProbes2$AF > 0.01), ]
snpProbes3 <- snpProbes3[which(snpProbes3$AF > 0.01), ]
snpProbes4 <- snpProbes4[which(snpProbes4$AF > 0.01), ]

dist <- cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist <- apply(dist, 1, min)
snpProbes4 <- snpProbes4[which(dist <= 10), ]

remove <- unique(c(tofilter$IlmnID, crosshyb[, 1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))
remove <- intersect(remove, rownames(res[[1]]))

if (length(remove) > 0) {
    for (i in 1:3) {
        res[[i]] <- res[[i]][-match(remove, rownames(res[[i]])), ]
    }
}

# add gene annotation
probeAnnot <- read.table(file.path(refPath, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
probeAnnot <- probeAnnot[match(rownames(res[[1]]), probeAnnot$probeID), ]
probeAnnot[["chrm"]] <- gsub("chr", "", probeAnnot[["chrm"]], fixed = TRUE)
probeAnnot$start <- probeAnnot$start+1

for(i in 1:3){
    res[[i]]<- cbind(res[[i]], probeAnnot[, c("chrm", "start", "GeneNames", "GeneClasses", "CGI", "CGIPosition")])
    res[[i]]<- res[[i]][!res[[i]][,"chrm"] %in% c("Y", "*"),]
}

# load GWAS regions

gwasClumped<-read.csv(gwasFile, stringsAsFactors = FALSE)
#convert chr X to 23
gwasClumped$CHR[which(gwasClumped$CHR == 23)]<-"X"


#----------------------------------------------------------------------#
# IDENTIFY PROBES IN REGIONS
#----------------------------------------------------------------------#

gwasRegions<-GRanges(seqnames = gwasClumped$CHR, strand = "*", 
        ranges = IRanges(start = gwasClumped$range.left, end = gwasClumped$range.right))
ewasRegions<-GRanges(seqnames = res[[1]]$chrm, strand = "*", 
        ranges = IRanges(start = res[[1]]$start, end = res[[1]]$start))
mcols(ewasRegions)<-cbind(res[[1]][, c("NullModel_SCZ_P", "NullModel_SCZ_coeff")],
res[[2]][, c("NullModel_SCZ_P", "NullModel_SCZ_coeff")],
res[[3]][, c("NullModel_SCZ_P", "NullModel_SCZ_coeff")], rownames(res[[1]]))
colnames(mcols(ewasRegions))<-c(outer(c("P", "Coeff"), cellTypes, paste, sep = ":"), "ProbeID")


interSect<-findOverlaps(ewasRegions, gwasRegions)
indexAnyRegion<-unique(queryHits(interSect))

#----------------------------------------------------------------------#
# TEST ENRICHMENT
#----------------------------------------------------------------------#

# Treat as a "pathway": evidence of smaller p-values?
# Any discovery DMPs in these regions?
# does being in a region predict being a DMP? 

thres<-1e-6
mwP <- rep(NA, 3)
nDMPs <- rep(NA, 3)
regressionP <- rep(NA, 3)
region_ind <- rep(0, length(ewasRegions))
region_ind[indexAnyRegion] <- 1

for(i in seq(1, 6, 2)){
    mwP[(i+1)/2] <- wilcox.test(mcols(ewasRegions)[indexAnyRegion,i], mcols(ewasRegions)[-indexAnyRegion,i])$p.value
    nDMPs[(i+1)/2] <- sum(mcols(ewasRegions)[indexAnyRegion,i] < thres)
    DMP_ind <- as.numeric(mcols(ewasRegions)[,i] < thres)
    if(sum(DMP_ind) > 0){
        model<- glm(DMP_ind ~ region_ind)
        regressionP[(i+1)/2] <- summary(model)$coefficients["region_ind",4]

        # save DMPs that overlap a GWAS region
        write.csv(res[[(i+1)/2]][DMP_ind & region_ind,], 
            file.path(resPath, "Tables", paste0("Discovery", cellTypes[(i+1)/2],"DMPsinSCZGWASRegions.csv")))
    }
}

outSum <- rbind(mwP, regressionP, nDMPs)
rownames(outSum) <- c("MannWhitneyP", "RegressionP", "nDMPs")
colnames(outSum)<- cellTypes

write.csv(outSum,
            file.path(resPath, "Tables", paste0("SummaryStatisticsEnrichmentofDMPsinSCZGWASRegionsAll.csv")))

#----------------------------------------------------------------------#
# CALCULATE COMBINED P-VALUES FOR EACH REGION
#----------------------------------------------------------------------#

probeIndex<-queryHits(interSect)
regionIndex<-subjectHits(interSect)

out<-matrix(data = NA, nrow = nrow(gwasClumped), ncol = 19)
colnames(out)<-c("nProbe", c(t(outer(cellTypes, c("nDMPs", "minCaseConP", "FishersCaseConP", "MeanCor", "MaxCor", "BrownsCaseConP"), 
    paste0))))

for(i in unique(regionIndex)){
	res.sub<-mcols(ewasRegions)[probeIndex[which(regionIndex == i)],]
	nProbes<-nrow(res.sub)
	out[i,1]<-nProbes
    outStat<-NULL
	if(nProbes > 1){
        for(j in seq(1, 6, 2)){
           pval<-res.sub[,j]
           # pvalue sum stats
    	    outStat<-c(outStat, sum(pval < thres),
	            min(pval))
            # calc fisher's combined p
		    fisher<-(sum(-log(pval))*2)
		    outStat<-c(outStat, 1-pchisq(fisher, 2*length(pval)))
			#calc brown's combined p
			betas.sub<-celltypeNormbeta[res.sub$ProbeID,which(QCmetrics$Cell.type == cellTypes[(j+1)/2])]
			corMat<-cor(t(betas.sub))
			cor.vector<-c(unique(as.vector(corMat)))[-1]
			covar<-calcCovMatrix(cor.vector)
			brownsp<-brownsP(covar, pval)
			outStat<-c(outStat, mean(abs(cor.vector)), max(cor.vector), brownsp[1])
		}
	} else {
        outStat<-rep(NA, 18)
    }
    out[i,-1]<-outStat
} 

out<-cbind(gwasClumped[,c("SNP","CHR","P","OR" ,"range.left","range.right")], out)
write.csv(out, file.path(resPath, "Tables", "CombinedPvaluesForSCZGWASRegionsTrubetskoyetal.csv"), 
    row.names = FALSE)



