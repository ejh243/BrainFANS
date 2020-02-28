## Written by Eilis
## takes .info files (1 per chromosome) and creates plots to summarise imputation


args<-commandArgs(trailingOnly = TRUE)

library(data.table)
dirName<-args[1]
refFile<-args[2]
colName<-args[3]

refPanel<-fread(refFile, data.table = FALSE)
if(!"id" %in% colnames(refPanel)){
	refPanel$id<-paste(refPanel$"#CHROM", refPanel$POS, refPanel$REF, refPanel$ALT, sep = ":")
}


compRefPanel<-matrix(data = NA, ncol = 2, nrow = 22)
colnames(compRefPanel)<-c("cor", "MAD")
nVar<-matrix(data = 0, nrow = 10, ncol = 10)
for(chr in 1:22){
	png(paste(dirName,"/ImputationQualityPlots_chr", chr, ".png", sep = ""), width = 12, height = 8, units = "in", res = 200)
	par(mfrow = c(2,2))
	
	imputScores<-read.table(gzfile(paste(dirName,"/chr", chr, ".info.gz", sep = "")), header = TRUE, na.strings = "-")
	index<-match(imputScores$SNP,refPanel$id)

	plot(density(imputScores$Rsq, na.rm = TRUE), main = paste("chr", chr, sep = ""), xlab = "Rsq")

	## compare freq against refPanel
	plot(refPanel[index,colName],imputScores$ALT_Frq,xlab = "RefPanel MAF", ylab = "Imputed sample MAF", pch = 16)
	## calculate median absolute deviation
	mtext(side = 3, line = 0.5, adj = 1, paste("MAD =", signif(median(abs(refPanel[index,colName] -imputScores$ALT_Frq), na.rm = TRUE), 3), sep = ""))
	
	compRefPanel[chr,1]<-cor(refPanel[index,colName],imputScores$ALT_Frq, use = "p")
	compRefPanel[chr,2]<-median(abs(refPanel[index,colName] -imputScores$ALT_Frq), na.rm = TRUE)	

	## summarise INFO score distribution by MAF	
	tabRsqMAF<-table(cut(imputScores$Rsq, seq(0,1,0.1)), cut(imputScores$MAF, seq(0,0.5,0.05)))
	barplot(tabRsqMAF,ylab = "N variants", legend = rownames(tabRsqMAF), xlab = "MAF", col = colorRampPalette(c("white", "navy"))(10), main = paste("chr", chr, sep = ""))
	tabRsqMAFPer<-t(t(tabRsqMAF)/colSums(tabRsqMAF)*100)
	barplot(tabRsqMAFPer,ylab = "%", xlab = "MAF", col = colorRampPalette(c("white", "navy"))(10))

	nVar<-nVar+tabRsqMAF
	

	dev.off()
}

write.csv(nVar, paste(dirName,"/nVariantsSummary.csv", sep = ""))
write.csv(compRefPanel, paste(dirName,"/CompRefPanelSummary.csv", sep = ""))

