
##---------------------------------------------------------------------#
##
## Title: Calculate data quality metrics from raw DNAm data 
##
## Purpose of script: From GDS file generate summary metrics for stages 1 & 2 of quality control filtering
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# project folder is provided on command line
# path to folder where cortical clock coefficients are saved is provided on command line
# excludes samples with really low intensities values (< 500) at beginning
# requires gdsfile to already be generated
# assumes matched genotype data, if it exists is in the 0_metadata folder named epicSNPs.raw

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refDir <- args[2]

gdsFile <-paste0(dataDir, "/2_gds/raw.gds")
qcData <-paste0(dataDir, "/2_gds/QCmetrics/QCmetrics.rdata")
genoFile <- paste0(dataDir, "/0_metadata/epicSNPs.raw")
configFile <- paste0(dataDir, "/config.r")

source(configFile)

arrayType <- toupper(arrayType)

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(e1071, warn.conflicts = FALSE, quietly = TRUE)
library(data.table, warn.conflicts = FALSE, quietly = TRUE)
library(bigmelon, warn.conflicts = FALSE, quietly = TRUE)
library(wateRmelon, warn.conflicts = FALSE, quietly = TRUE)


#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#
setwd(dataDir)

# load sample sheet
sampleSheet<-read.csv("0_metadata/sampleSheet.csv", na.strings = c("", "NA"), stringsAsFactors = FALSE)
# if no column Basename, creates from columns Chip.ID and Chip.Location
if(!"Basename" %in% colnames(sampleSheet)){
	sampleSheet$Basename<-paste(sampleSheet$Chip.ID, sampleSheet$Chip.Location, sep = "_")
}


gfile<-openfn.gds(gdsFile, readonly = FALSE, allow.fork = TRUE)
# ensure sample sheet is in same order as data
sampleSheet<-sampleSheet[match(colnames(gfile), sampleSheet$Basename),]

## see if any QC data already exists
if(file.exists(qcData)){
	load(qcData)
	## check contains all required samples
	if(nrow(QCmetrics) == nrow(sampleSheet)){
		print("QC file loaded")
	} else {
		QCmetrics<-sampleSheet
		print("QC file to be updated with new samples")
	}
} else{
	QCmetrics<-sampleSheet
	print("QC object initiated")
}

manifest <- cdegUtilities::readManifest(
	referenceDirectory = refDir,
	probeMatchingIndex = fData(gfile)[["Probe_ID"]],
	arrayType = arrayType 
)
if (!exists("manifest"))
	stop("Manifest file could not be loaded correctly")



#----------------------------------------------------------------------#
# CALCULATE QC METRICS
#----------------------------------------------------------------------#

# calculate median M & U intensity
if(!"M.median" %in% colnames(QCmetrics)){
	print("Calculating signal intensity statistics")

	m_intensities<-methylated(gfile)
	u_intensities<-unmethylated(gfile)
	M.median<-unlist(apply.gdsn(m_intensities, 2, median, na.rm = TRUE))
	U.median<-unlist(apply.gdsn(u_intensities, 2, median, na.rm = TRUE))

	intens.ratio<-M.median/U.median
	# exclude really poor intensity samples from beginning so rest of QC is not dominated by them
	intensPASS<-M.median > 500
	# exclude fully methylated control samples
	intensPASS[which(intens.ratio > 4)]<-FALSE
	QCmetrics<-cbind(QCmetrics,M.median, U.median, intens.ratio, intensPASS)

} else {
	intensPASS<-QCmetrics$intensPASS
}


# calculate bisulfite conversion statistic
if(!"bisulfCon" %in% colnames(QCmetrics)){	
	print("Calculating bisulfite conversion statistics")

	bisulfCon<-bscon(gfile)

	bisulfCon[which(intensPASS == FALSE)]<-NA
	QCmetrics<-cbind(QCmetrics,bisulfCon)
}


# PCA of control-probe intensities
if(!"PC1_cp" %in% colnames(QCmetrics)){	
	print("Calculating PCs of control probes")
	# exclude really poor intensity samples

	qc.meth<-QCmethylated(gfile)[,QCmetrics$intensPASS]
	qc.unmeth<-QCunmethylated(gfile)[,QCmetrics$intensPASS]


	# remove negative controls
	qc.meth<-qc.meth[grep("Negative", rownames(qc.meth), invert=TRUE),]
	qc.unmeth<-qc.unmeth[grep("Negative", rownames(qc.unmeth), invert=TRUE),]
	ctrl.all<-t(rbind(qc.meth, qc.unmeth))

	# exclude columns where all NAs
	ctrl.all<-ctrl.all[,which(colSums(is.na(ctrl.all)) < nrow(ctrl.all))]

	pca <- prcomp(na.omit(ctrl.all))
	ctrlprobes.scores = pca$x
	colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
	ctrl.pca<-pca$sdev^2/sum(pca$sdev^2)
	ctrlprobes.scores<-ctrlprobes.scores[match(QCmetrics$Basename, QCmetrics$Basename[QCmetrics$intensPASS]),]
	rownames(ctrlprobes.scores)<-QCmetrics$Basename	
	# only save PCs which explain > 1% of the variance
	QCmetrics<-cbind(QCmetrics,ctrlprobes.scores[,which(ctrl.pca > 0.01)])
}


# perform PCA on beta values
if(!"PC1_betas" %in% colnames(QCmetrics)){
	print("Calculating PCs of autosomal beta values")
	# filter to autosomal only
	if(arrayType == "V2" | arrayType == "450K"){
		auto.probes<-which(manifest$CHR != "chrX" & manifest$CHR != "chrY")
	} else {
		auto.probes<-which(fData(gfile)$chr != "chrX" & fData(gfile)$chr != "chrY")
	}

	pca <- prcomp(t(na.omit(betas(gfile)[,][auto.probes,QCmetrics$intensPASS])))
	betas.scores = pca$x
	colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
	betas.pca<-pca$sdev^2/sum(pca$sdev^2)
	betas.scores<-betas.scores[match(QCmetrics$Basename, QCmetrics$Basename[QCmetrics$intensPASS]),]
	rownames(betas.scores)<-QCmetrics$Basename	
	# only save PCs which explain > 1% of the variance
	QCmetrics<-cbind(QCmetrics,betas.scores[,which(betas.pca > 0.01)])
}

# Identify outlier samples 
# excluded as comparable to PC filtering already included
#if(!"iqr" %in% colnames(QCmetrics)){
#outlierDetect<- outlyx(rawbetas, plot = FALSE)
#QCmetrics<-cbind(QCmetrics,outlierDetect)
#}


# detection p value filtering for sample filtering
if(!"pFilter" %in% colnames(QCmetrics)){	
	print("Running pfilter at sample level")

	pFOut<-apply.gdsn(
		node = pvals(gfile),
		margin = 2,
		FUN = function(x, y, z) {
			(sum(x > y, na.rm = TRUE)) < ((sum(!is.na(x)) * z)/100)
		},
		as.is = "logical",
		y = 0.05,
		z = 1
	)

	pFOut[!QCmetrics$intensPASS]<-NA
	QCmetrics<-cbind(QCmetrics,"pFilter"= pFOut)
}


# detection p value and beacd count filtering for probe filtering
if(!exists("probeFilt")){
	# exclude really poor intensity, low bisulf conversion and failed pFilt samples
	goodsamps <- QCmetrics$intensPASS & QCmetrics$bisulfCon > thresBS & QCmetrics$pFilter

	print("Running pfilter at probe level")
	pFiltProbesPass<-apply.gdsn(node = pvals(gfile),
	 margin = 1, FUN = function(x, y, z, passed) {
				(sum(x[passed] > y, na.rm = TRUE)) < ((sum(!is.na(x)) * z)/100)
			}, as.is = "logical", y = 0.05, z = 1, passed = goodsamps)

	print("Calculating probe beadcounts")
	beadcounts <- apply(index.gdsn(gfile, 'NBeads')[,], MARGIN = 1,
                            FUN = function(x,passed){
                            x[x<3] <- NA
                           	length(which(is.na(x[passed])=="TRUE"))
                            }, passed = goodsamps)
	beadFiltPass <- beadcounts<((nrow(QCmetrics) * 1)/100)

	probeFilt <- cbind(pFiltProbesPass, beadFiltPass)	
}


# NOTE this currently averages beta values accross duplicated probes
# calc Horvaths epigenetic age
if(!"DNAmAge" %in% colnames(QCmetrics)){
	print("Calculating Horvath's pan tissue epigenetic age")
	if(arrayType == "V2"){
		DNAmAge<-agep(epicv2clean(betas(gfile)[]))
	} else {
		data(coef)
		DNAmAge<-agep(gfile, coef=coef)
	}
	colnames(DNAmAge)[1] <- "DNAmAge"
	DNAmAge[!QCmetrics$intensPASS,]<-NA
	QCmetrics<-cbind(QCmetrics,DNAmAge)
}  

# NOTE this currently averages beta values accross duplicated probes
# calc Cortical Clock Age
if ((!"CCDNAmAge" %in% colnames(QCmetrics)) && tolower(tissueType) == "brain") {
	print("Calculating Shireby's Cortical Clock epigenetic age")
	CC_coef<-read.csv(paste0(refDir, "/CortexClock/CorticalClockCoefficients.csv"), stringsAsFactors = FALSE)
	anti.trafo <- function(x,adult.age=20) { 
		ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) 
	}
	if(arrayType == "V2"){
		cc <- CC_coef[CC_coef$probe %in% row.names(epicv2clean(betas(gfile)[])),]
		CCDNAmAge<- anti.trafo(as.numeric(CC_coef[1,2] + t(epicv2clean(betas(gfile)[])[row.names(epicv2clean(betas(gfile)[])) %in% CC_coef[-1,1],]) %*% cc[,2]))
	} else {
		CCDNAmAge<- anti.trafo(as.numeric(CC_coef[1,2] + t(betas(gfile)[][CC_coef[-1,1],]) %*% CC_coef[-1,2]))
	}

	QCmetrics<-cbind(QCmetrics,CCDNAmAge)

}




# check effect of normalisation
if(!"rmsd" %in% colnames(QCmetrics)){
	print("Calculating effect of normalisation")

	normbeta <- adjustedDasen(
		onetwo = manifest$designType,
		chr = manifest$CHR,
		mns = read.gdsn(methylated(gfile)),
		uns = read.gdsn(unmethylated(gfile)))
	add.gdsn(gfile, "normbeta", normbeta, replace = TRUE)
	
	qualDat<-qual(betas(gfile)[,], normbeta)
	qualDat[which(intensPASS == FALSE),]<-NA
	QCmetrics<-cbind(QCmetrics,qualDat)
}

# count number of missing values
if(!"nNAsPer" %in% colnames(QCmetrics)){
	print("Counting the number of missing beta values per sample")
	nNAs<-colSums(is.na(betas(gfile)[,]))
	nNAsPer<-nNAs/nrow(betas(gfile)[,])*100
	QCmetrics<-cbind(QCmetrics,nNAs, nNAsPer)
}

#----------------------------------------------------------------------#
# PREDICT SEX
#----------------------------------------------------------------------#

# NOTE threshold for M prediction not valid for epicV2 data
if(!"predSex" %in% colnames(QCmetrics)){	
	print("Performing sex prediction from sex chromosome profiles")	
	if(arrayType == "V2" | arrayType == "450K"){
		x.probes<-which(manifest$CHR == "chrX")
		y.probes<-which(manifest$CHR == "chrY")
	} else {
		x.probes<-which(fData(gfile)$chr == "chrX")
		y.probes<-which(fData(gfile)$chr == "chrY")
	}
	ints.auto<-methylated(gfile)[c(x.probes, y.probes),]+unmethylated(gfile)[c(x.probes, y.probes),]
	ints.X<-methylated(gfile)[x.probes,]+unmethylated(gfile)[x.probes,]
	ints.Y<-methylated(gfile)[y.probes,]+unmethylated(gfile)[y.probes,]

	x.cp<-colMeans(ints.X, na.rm = TRUE)/colMeans(ints.auto, na.rm = TRUE)
	y.cp<-colMeans(ints.Y, na.rm = TRUE)/colMeans(ints.auto, na.rm = TRUE)

	# test for evidence of two groups (i.e. catch for if dataset only has one sex)
	# result is identical for both chromosomes
	# needs at least ~15 samples for this to work
	library(diptest, warn.conflicts = FALSE, quietly = TRUE)
	bimodP.y<-dip.test(y.cp)$p.value

	if(bimodP.y < 0.05){
		print("Intensities on X & Y chromosomes are multimodal. Sufficient evidence of more than one sex in data, proceeding to fit two distributions to the data to make predictions")

		library(mixtools, warn.conflicts = FALSE, quietly = TRUE)
		mixmdl.x = normalmixEM(x.cp, mu = xMus, sigma = xSigmas)
		mixmdl.y = normalmixEM(y.cp, mu = yMus, sigma = ySigmas)

		posteriorProbs<-cbind(mixmdl.x$posterior, mixmdl.y$posterior)
		colnames(posteriorProbs)<-c("PP.M.X", "PP.F.X", "PP.F.Y", "PP.M.Y")
		# base prediction on y chromosome posterior probability
		predSex.y<-rep(NA, length(y.cp))
		predSex.y[which(posteriorProbs[,"PP.M.Y"] > 0.5 & intensPASS == TRUE)]<-"M"
		predSex.y[which(posteriorProbs[,"PP.F.Y"] > 0.5 & intensPASS == TRUE)]<-"F"

		# base prediction on x chromosome posterior probability
		predSex.x<-rep(NA, length(x.cp))
		predSex.x[which(posteriorProbs[,"PP.M.X"] > 0.5 & intensPASS == TRUE)]<-"M"
		predSex.x[which(posteriorProbs[,"PP.F.X"] > 0.5 & intensPASS == TRUE)]<-"F"

		# check for consistent prediction
		predSex<-rep(NA, length(x.cp))
		predSex[which(predSex.x == predSex.y)]<-predSex.x[which(predSex.x == predSex.y)]
		QCmetrics<-cbind(QCmetrics,x.cp,y.cp,posteriorProbs,predSex.x, predSex.y, predSex)
	} else{
		print("Intensities on X & Y chromosomes are unimodal. This would suggest just one sex in the data, a sample sample size or very few of one sex. Bypassing prediction step")
		predSex<-rep(NA, length(x.cp))
		predSex.x<-rep(NA, length(x.cp))
		predSex.y<-rep(NA, length(x.cp))
		QCmetrics<-cbind(QCmetrics,x.cp,y.cp,predSex.x, predSex.y, predSex)

	}
}

#----------------------------------------------------------------------#
# CORRELATION ACROSS ARRAY SNPS
#----------------------------------------------------------------------#

# check duplicate samples using SNPs on array
if(!exists("snpCor")){
	print("Calculating pairwise correlations across SNP probes")	
	rsbetas<-betas(gfile)[,][grep("rs", rownames(betas(gfile)[,])),]
	snpCor<-cor(rsbetas, use = "pairwise.complete.obs")
}

#----------------------------------------------------------------------#
# COMPARE TO EXTERNAL SNP DATA
#----------------------------------------------------------------------#

indexIID <- grep("^Genotype_IID$", names(QCmetrics), ignore.case=TRUE)
if(!"genoCheck"%in% colnames(QCmetrics) & file.exists(genoFile)){
	print("Comparing against matched genotype data")
	geno<-read.table(genoFile, stringsAsFactors = FALSE, header = TRUE)
	geno.all<-geno

	if (exists("indexIID") & length(indexIID) != 1){
		message("Warning: Genotype_IID column is missing, unable to compare external SNP data.")
	}else{
		geno<-geno[match(QCmetrics$Genotype_IID, geno$IID),]
		rsIDs<-gsub("_.", "", colnames(geno)[-c(1:6)])

		if(arrayType == "V2"){
			betas.rs<-epicv2clean(betas(gfile)[,])[rsIDs,]
		} else {
			betas.rs<-betas(gfile)[,][rsIDs,]
		}


		# first check direction of minor alleles
		cors<-vector(length = length(rsIDs))
		for(i in 1:length(rsIDs)){
			cors[i]<-cor(geno[,i+6], betas.rs[i,], use = "pairwise.complete.obs")
		}
		# change minor allele in genotype data if negative correlation
		for(each in which(cors < 0)){
			geno[,each+6]<-(2-geno[,each+6])
			geno.all[,each+6]<-(2-geno.all[,each+6])
		}
		geno.mat<-as.matrix(geno[,-c(1:6)])
		geno.all.mat<-as.matrix(geno.all[,-c(1:6)])
		rownames(geno.all.mat)<-geno.all$IID

		genoCheck<-rep(NA, nrow(QCmetrics))
		for(i in 1:ncol(betas.rs)){
			if(!is.na(geno[i,1]) & QCmetrics$intensPASS[i] == TRUE){
				genoCheck[i]<-cor(geno.mat[i,], betas.rs[,i], use = "pairwise.complete.obs")
			}
		}

		# if any incongruent perform search for best using all geno data
		# first though check if any geno combinations present multiple times in this cohort:
		# count how many individuals with each geno combination in sample
		indGenoCombo<-apply(geno.all.mat, 1, paste, collapse = ";")
		tabGeneticInd<-table(indGenoCombo)
		#table(tabGeneticInd)

		# pull out list of samples which identical genotypes across these variants
		dupCombos<-names(tabGeneticInd[which(tabGeneticInd > 1)])
		if(length(dupCombos) > 0){
			dupIDs<-NULL
			for(each in dupCombos){
				dupIDs<-c(dupIDs, paste(rownames(geno.all.mat)[which(indGenoCombo == each)], collapse = ";"))
			}
			write.csv(dupIDs, paste0(dataDir, "/2_gds/QCmetrics/IndividualsWithIdenticalGenotypeCombinationsInComparisionWithSNPData.csv"))
		}


		genoMatch<-rep(NA, nrow(QCmetrics))
		genoMatchVal<-rep(NA, nrow(QCmetrics))
		for(i in 1:ncol(betas.rs)){
			if(QCmetrics$intensPASS[i] == TRUE){
				corVals<-rep(NA, nrow(geno.all.mat))
				for(j in 1:nrow(geno.all.mat)){
					corVals[j]<-cor(geno.all.mat[j,], betas.rs[,i], use = "pairwise.complete.obs")
				}
			}
			if(max(corVals, na.rm = TRUE) > 0.9){ 
				# as possible to match multipe, save all
				genoMatch[i]<-paste(rownames(geno.all.mat)[which(corVals > 0.9)], collapse = ";")
				genoMatchVal[i]<-paste(corVals[which(corVals > 0.9)], collapse = ";")
			}
		}
		QCmetrics<-cbind(QCmetrics,genoCheck, genoMatch, genoMatchVal)
	}

}



## perform sample filtering

QCSum<-cbind(QCmetrics$bisulfCon > thresBS,(QCmetrics$M.median > intenThres & QCmetrics$U.median > intenThres), QCmetrics$rmsd < nvThres, QCmetrics$nNAsPer < perMiss, QCmetrics$pFilter, QCmetrics$intens.ratio < 4)
colnames(QCSum)<-c(paste0("BSConversion>",thresBS),  paste0("intens>", intenThres), paste0("normviolence<", nvThres), paste0("%MissingSites<", perMiss), "DetectionPvalues", "MethylatedControl")

QCSum<-cbind(QCSum, rowSums(QCSum, na.rm = TRUE) == rowSums(!is.na(QCSum)))
colnames(QCSum)[ncol(QCSum)]<-"passQCS1"

if(sexCheck){
	sexPass <- as.character(QCmetrics$predSex) == as.character(QCmetrics$Sex)
	sexPass[is.na(QCmetrics$predSex)]<-FALSE
	## if unable to predict sex, change to fail
	QCSum<-cbind(QCSum, sexPass)
}else{
	QCSum<-cbind(QCSum, NA)
}
colnames(QCSum)[ncol(QCSum)]<-"sexCheck"

if(snpCheck & "genoCheck"%in% colnames(QCmetrics)){
	QCSum<-cbind(QCSum, QCmetrics$genoCheck > 0.8)
}else{
	QCSum<-cbind(QCSum, NA)
}
colnames(QCSum)[ncol(QCSum)]<-"snpCheck"

# if neither sex check or snp check performed stage two pass status is copied from stage one 
cols<-c("passQCS1", "sexCheck", "snpCheck")
cols<-cols[cols %in% colnames(QCSum)]
if(length(cols) > 1){
	QCSum<-cbind(QCSum, rowSums(QCSum[,cols], na.rm = TRUE) == rowSums(!is.na(QCSum[,cols])))
} else {
	QCSum<-cbind(QCSum, QCSum[,"passQCS1"])
}
colnames(QCSum)[ncol(QCSum)]<-"passQCS2"
rownames(QCSum)<-QCmetrics$Basename


sampleColsToKeep<-c("Basename", "Sample_ID", "Individual_ID", "Cell_Type")
sampleColsToKeep<-sampleColsToKeep[sampleColsToKeep %in% colnames(QCmetrics)]
write.csv(cbind(QCmetrics[,sampleColsToKeep], QCSum),  paste0(dataDir, "/2_gds/QCmetrics/PassQCStatusAllSamples.csv"))


#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

closefn.gds(gfile)

write.csv(QCmetrics, paste0(dataDir, "/2_gds/QCmetrics/QCMetricsPostSampleCheck.csv"))

# save QC metrics and SNP correlations to generate QC report
if(file.exists(genoFile) & exists("indexIID") & length(indexIID) == 1 ){
	save(QCmetrics, probeFilt, snpCor, betas.pca, ctrl.pca, pFOut, geno.mat, betas.rs, rsbetas, file = qcData)
} else {
	save(QCmetrics, probeFilt, snpCor, betas.pca, ctrl.pca, pFOut, rsbetas, file = qcData)
}

