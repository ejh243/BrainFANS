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

# calculate median M & U intensity
if(!"M.median" %in% colnames(QCmetrics)){
	
	m_intensities<-methylated(gfile)
	u_intensities<-unmethylated(gfile)
	M.median<-apply.gdsn(m_intensities, 2, median)
	U.median<-apply.gdsn(u_intensities, 2, median)
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
	QCmetrics<-cbind(QCmetrics,ctrlprobes.scores)
}

# Identiify outlier samples
rawbetas<-betas(gfile)
outlyx<- outlyx(rawbetas, plot = FALSE)
QCmetrics<-cbind(QCmetrics,outlyx)

# detection p value filtering
pfilter.gds(gfile)


# calc Horvaths epigenetic age
data(coef)
DNAmAge<-agep(gfile, coef=coef)
QCmetrics<-cbind(QCmetrics,DNAmAge)

# check sex
predictSex(gfile)

# check effect of normalisation


## need to close gds file in order to open in another R session
closefn.gds(gfile)


save(QCmetrics, file = qcData)