## run this script to calculate QC metrics
## requires gdsfile to already be generated and sample sheet to be loaded

setwd(dataDir)

gfile<-openfn.gds(gdsFile, readonly = FALSE)

## see if any QC data already exists
if(file.exists(qcData)){
	load(qcData)
} else{
	QCmetrics<-sampleSheet
}

## extract a few useful matrices
rawbetas<-betas(gfile)[,]

# calculate median M & U intensity
if(!"M.median" %in% colnames(QCmetrics)){
	
	m_intensities<-methylated(gfile)
	u_intensities<-unmethylated(gfile)
	M.median<-unlist(apply.gdsn(m_intensities, 2, median))
	U.median<-unlist(apply.gdsn(u_intensities, 2, median))
	QCmetrics<-cbind(QCmetrics,M.median, U.median)
	
}

# calculate bisulfite conversion statistic
if(!"bisulfCon" %in% colnames(QCmetrics)){	
	bisulfCon<-bscon(gfile)
	QCmetrics<-cbind(QCmetrics,bisulfCon)
}

## PCA of control-probe intensities
if(!"PC1_cp" %in% colnames(QCmetrics)){	
	qc.meth<-QCmethylated(gfile)
	qc.unmeth<-QCunmethylated(gfile)
	# remove negative controls
	qc.meth<-qc.meth[grep("Negative", rownames(qc.meth), invert=TRUE),]
	qc.unmeth<-qc.unmeth[grep("Negative", rownames(qc.unmeth), invert=TRUE),]
	ctrl.all<-t(rbind(qc.meth, qc.unmeth))

	pca <- prcomp(na.omit(ctrl.all))
	ctrlprobes.scores = pca$x
	colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
	## only store the first 20
	QCmetrics<-cbind(QCmetrics,ctrlprobes.scores[,1:20])
}

## perform PCA on beta values
if(!"PC1_betas" %in% colnames(QCmetrics)){	

	pca <- prcomp(na.omit(rawbetas))
	betas.scores = pca$x
	colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
	## only store the first 20
	QCmetrics<-cbind(QCmetrics,betas.scores[,1:20])
}

## Identiify outlier samples
if(!"iqr" %in% colnames(QCmetrics)){
	outlierDetect<- outlyx(rawbetas, plot = FALSE)
	QCmetrics<-cbind(QCmetrics,outlierDetect)
}

## detection p value filtering
pfilter.gds(gfile)



## calc Horvaths epigenetic age
if(!"DNAmAge" %in% colnames(QCmetrics)){	
	data(coef)
	DNAmAge<-agep(gfile, coef=coef)
	QCmetrics<-cbind(QCmetrics,DNAmAge)
}

## check sex
# idenitfy X chromosome probes
if(!"predSex" %in% colnames(QCmetrics)){	
	x.probes<-which(fData(gfile)$chr == "chrX")
	predSex<-predictSex(gfile, x.probes)
	QCmetrics<-cbind(QCmetrics,predSex)
}

## check duplicate samples using SNPs on array
if(!exists("snpCor")){
	rsbetas<-rawbetas[grep("rs", rownames(rawbetas)),]
	snpCor<-cor(rsbetas, use = "pairwise.complete.obs")
}

# check effect of normalisation
if(!"rmsd" %in% colnames(QCmetrics)){
	dasen(gfile, node="normbeta")
	normbetas<-index.gdsn(gfile, "normbeta")[,]
	qualDat<-qual(rawbetas, normbetas)
	QCmetrics<-cbind(QCmetrics,qualDat)
}

## need to close gds file in order to open in another R session
closefn.gds(gfile)

# save QC metrics and SNP correlations to generate QC report
save(QCmetrics, snpCor, file = qcData)