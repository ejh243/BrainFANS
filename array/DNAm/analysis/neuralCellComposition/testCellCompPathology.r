##---------------------------------------------------------------------#
##
## Title: Test estimated cell composition against AD neuropathology
##
## Purpose of script: Test estimated cell composition against AD neuropathology
##
## Author: Eilis Hannon
##
## Date Created: 13/12/2022
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# loads cell composition predicions from rdata objects


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
## set plotting colours
library(paletteer)
ctCols <- paletteer_d("ggsci::category10_d3")
names(ctCols)<-c("DoubleNeg", "NeuNPos", "NeuNNeg", "Sox10Pos", "IRF8Pos", "TripleNeg", "SATB2Neg","SATB2Pos", 
"SOX6Neg", "SOX6Pos")

args<-commandArgs(trailingOnly = TRUE)
bulkPath <- args[1]
plotPath <- args[2]

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(meta)

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

load(bulkPath)
pos <- position_dodge(0.9)

#----------------------------------------------------------------------#
# SUMMARISE DATASET
#----------------------------------------------------------------------#

sumStat<-cbind(table(cellCompAll[[1]]$Dataset), 
aggregate(Age ~ Dataset, cellCompAll[[1]], summary),
aggregate(Braak ~ Dataset, cellCompAll[[1]], summary))

write.csv(sumStat, file.path(plotPath, "SummariseDatasets.csv"))


#----------------------------------------------------------------------#
# QC CELL COMPOSITION ESTIMATES
#----------------------------------------------------------------------#
fig1a<-list()
fig1b<-list()
fig2a<-list()
fig2b<-list()
fig3a<-list()
fig3b<-list()

resAD<-NULL
resBraak<-NULL
resMeta<-NULL
for(i in 1:length(cellCompAll)){
	## define AD status
	cellCompAll[[i]]$AD<-rep(NA, nrow(cellCompAll[[i]]))
	cellCompAll[[i]]$AD[which(cellCompAll[[i]][,"Braak"] <3)]<-0
	cellCompAll[[i]]$AD[which(cellCompAll[[i]][,"Braak"] >4)]<-1
	cellCompAll[[i]]$AD <- as.factor(cellCompAll[[i]]$AD)
	
	# Violin plot of CETYGO
	fig1a[[i]]<-ggplot(subset(cellCompAll[[i]], Method == "ANOVA"), aes(x=Dataset, y=CETYGO, fill = Dataset))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "CETYGO")
					   
	fig1b[[i]]<-ggplot(subset(cellCompAll[[i]], Method == "IDOL"), aes(x=Dataset, y=CETYGO, fill = Dataset))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "CETYGO")
					   
	# CETYGO against Braak
	fig2a[[i]]<-ggplot(subset(cellCompAll[[i]], Method == "ANOVA" & !is.na(Braak)), aes(x=Dataset, y=CETYGO, fill = factor(Braak)))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "CETYGO")
		  
	fig2b[[i]]<-ggplot(subset(cellCompAll[[i]], Method == "IDOL" & !is.na(Braak)), aes(x=Dataset, y=CETYGO, fill = factor(Braak)))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "CETYGO")
		  

	# CETYGO against AD
	fig3a[[i]]<-ggplot(subset(cellCompAll[[i]], Method == "ANOVA" & !is.na(AD)), aes(x=Dataset, y=CETYGO, fill = AD))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "CETYGO")

	# CETYGO against AD
	fig3b[[i]]<-ggplot(subset(cellCompAll[[i]], Method == "IDOL" & !is.na(AD)), aes(x=Dataset, y=CETYGO, fill = AD))  +
		  geom_violin(position = pos, scale = 'width')  +
		  stat_summary(fun = "mean", 
					   geom = "point", 
					   position = pos, col = "white") + theme(legend.position = "none") + labs(x = "", y = "CETYGO")

		
	# Test CETYGO against AD & Braak
	for(method in c("ANOVA", "IDOL")){
		if(nrow(subset(cellCompAll[[i]], Method == method)) > 0){
		
			te<-NULL
			seTE<-NULL

			model<-lm(CETYGO ~ AD + Age + Sex + Dataset + AD*Dataset, data = subset(cellCompAll[[i]], Method == method))
			resAD<-rbind(resAD, c("Samples" = "All", "Method" = method, "Panel" = i, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))
				
			model<-lm(CETYGO ~ AD + Age + Sex , data = subset(cellCompAll[[i]], Method == method & Dataset == "BDR"))
			resAD<-rbind(resAD, c("Samples" = "BDR", "Method" = method, "Panel" = i, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))		
			te<-c(te, summary(model)$coefficients[c("AD1"),c(1)])
			seTE<-c(seTE, summary(model)$coefficients[c("AD1"),c(2)])				
			
			model<-lm(CETYGO ~ AD + Age + Sex , data = subset(cellCompAll[[i]], Method == method & Dataset == "LBB1"))
			resAD<-rbind(resAD, c("Samples" = "LBB1", "Method" = method, "Panel" = i, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))		
			te<-c(te, summary(model)$coefficients[c("AD1"),c(1)])
			seTE<-c(seTE, summary(model)$coefficients[c("AD1"),c(2)])
			
			model<-lm(CETYGO ~ AD + Age + Sex , data = subset(cellCompAll[[i]], Method == method & Dataset == "MS"))
			resAD<-rbind(resAD, c("Samples" = "MS", "Method" = method,"Panel" = i, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))		
			te<-c(te, summary(model)$coefficients[c("AD1"),c(1)])
			seTE<-c(seTE, summary(model)$coefficients[c("AD1"),c(2)])

			## meta analysis
			outMeta<-metagen(te, seTE)
			resMeta<-rbind(resMeta, c("AD", method, i,outMeta$TE.fixed,outMeta$seTE.fixed,outMeta$pval.fixed, outMeta$TE.random, outMeta$seTE.random,outMeta$pval.random, outMeta$tau, outMeta$I2, outMeta$Q,1-pchisq(outMeta$Q, outMeta$df.Q)))

			te<-NULL
			seTE<-NULL
			
			model<-lm(CETYGO ~ Braak + Age + Sex + Dataset + Braak*Dataset, data = subset(cellCompAll[[i]], Method == method))
			resBraak<-rbind(resBraak, c("Samples" = "All", "Method" = method, "Panel" = i, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))
				
			model<-lm(CETYGO ~ Braak + Age + Sex , data = subset(cellCompAll[[i]], Method == method & Dataset == "BDR"))
			resBraak<-rbind(resBraak, c("Samples" = "BDR", "Method" = method, "Panel" = i, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))		
			te<-c(te, summary(model)$coefficients[c("Braak"),c(1)])
			seTE<-c(seTE, summary(model)$coefficients[c("Braak"),c(2)])
					
			model<-lm(CETYGO ~ Braak + Age + Sex , data = subset(cellCompAll[[i]], Method == method & Dataset == "LBB1"))
			resBraak<-rbind(resBraak, c("Samples" = "LBB1", "Method" = method, "Panel" = i, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))	
			te<-c(te, summary(model)$coefficients[c("Braak"),c(1)])
			seTE<-c(seTE, summary(model)$coefficients[c("Braak"),c(2)])			
						
			model<-lm(CETYGO ~ Braak + Age + Sex , data = subset(cellCompAll[[i]], Method == method & Dataset == "MS"))
			resBraak<-rbind(resBraak, c("Samples" = "MS", "Method" = method,"Panel" = i, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))	
			te<-c(te, summary(model)$coefficients[c("Braak"),c(1)])
			seTE<-c(seTE, summary(model)$coefficients[c("Braak"),c(2)])
			
			outMeta<-metagen(te, seTE)
			resMeta<-rbind(resMeta, c("Braak", method, i,outMeta$TE.fixed,outMeta$seTE.fixed,outMeta$pval.fixed, outMeta$TE.random, outMeta$seTE.random,outMeta$pval.random, outMeta$tau, outMeta$I2, outMeta$Q,1-pchisq(outMeta$Q, outMeta$df.Q)))
			
			}			
	}
}

colnames(resAD)[-c(1:3)]<-c("ADCoeff", "ADSE", "ADP", "AgeCoeff", "AgeSE", "AgeP", "SexCoeff", "SexSE", "SexP")
colnames(resBraak)[-c(1:3)]<-c("BraakCoeff", "BraakSE", "BraakP", "AgeCoeff", "AgeSE", "AgeP", "SexCoeff", "SexSE", "SexP")
colnames(resMeta)<-c("Phenotype", "Method", "Panel", "FixedCoeff", "FixedSE", "FixedP", "RandomCoeff", "RandomSE", "RandomP", "tau", "I2", "Q", "HetP")

write.csv(resAD, file.path(plotPath, "CETYGOAgainstADStatus.csv"))
write.csv(resBraak, file.path(plotPath, "CETYGOAgainstBraakStatus.csv"))
write.csv(resMeta, file.path(plotPath, "CETYGOAgainstADStatusMetaAnalysis.csv"))


ggarrange(plotlist=fig1a, nrow = 2, ncol = 4)
ggsave(filename = file.path(dirname(bulkPath), "plots", "AD", "ViolinPlotCETYGOByDatasetANOVA.pdf"),  units = "in", width = 18, height = 8)

ggarrange(plotlist=fig2a, nrow = 2, ncol = 4)
ggsave(filename = file.path(dirname(bulkPath), "plots", "AD", "ViolinPlotCETYGOByBraakANOVA.pdf"),  units = "in", width = 18, height = 8)
	
ggarrange(plotlist=fig3a, nrow = 2, ncol = 4)
ggsave(filename = file.path(dirname(bulkPath), "plots", "AD", "ViolinPlotCETYGOByADStatusANOVA.pdf"),  units = "in", width = 18, height = 8)


ggarrange(plotlist=fig1b, nrow = 2, ncol = 4)
ggsave(filename = file.path(dirname(bulkPath), "plots", "AD", "ViolinPlotCETYGOByDatasetIDOL.pdf"),  units = "in", width = 18, height = 8)

ggarrange(plotlist=fig2b, nrow = 2, ncol = 4)
ggsave(filename = file.path(dirname(bulkPath), "plots", "AD", "ViolinPlotCETYGOByBraakIDOL.pdf"),  units = "in", width = 18, height = 8)
	
ggarrange(plotlist=fig3b, nrow = 2, ncol = 4)
ggsave(filename = file.path(dirname(bulkPath), "plots", "AD", "ViolinPlotCETYGOByADStatusIDOL.pdf"),  units = "in", width = 18, height = 8)

#----------------------------------------------------------------------#
# TEST AGAINST BRAAK AND AD
#----------------------------------------------------------------------#

## extract best prediction for each cell type

predCCBest<-cbind(subset(cellCompAll[[1]], Method == "ANOVA")[,c("DoubleNeg", "NeuNPos")], subset(cellCompAll[[5]], Method == "IDOL")[,c("IRF8Pos","TripleNeg")], subset(cellCompAll[[3]], Method == "ANOVA")[,c("SATB2Neg", "SATB2Pos")], subset(cellCompAll[[6]], Method == "ANOVA")[,c("NEUNNeg", "SOX6Pos", "SOX6Neg")], subset(cellCompAll[[4]], Method == "IDOL")[,"Sox10Pos"])
colnames(predCCBest)[10]<-"Sox10Pos"
predCCBest<-cbind(predCCBest, cellCompAll[[1]][,c("Dataset", "Age", "Sex", "Braak")])

predCCBest$Braak_factor <- as.factor(predCCBest$Braak)

resOut<-NULL
resMeta<-NULL
fig1a<-list()

for(j in 1:10){
	ct<-colnames(predCCBest)[j]

	te<-NULL
	seTE<-NULL

	model<-lm(get(ct) ~ Braak + Age + Sex + Dataset + Braak*Dataset, data = predCCBest)
	resOut<-rbind(resOut, c("Samples" = "All","CellType" = ct, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))
	
	model<-lm(get(ct) ~ Braak + Age + Sex , data = subset(predCCBest,Dataset == "BDR"))
	resOut<-rbind(resOut, c("Samples" = "BDR", "CellType" = ct, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))		
	te<-c(te, summary(model)$coefficients[c("Braak"),c(1)])
	seTE<-c(seTE, summary(model)$coefficients[c("Braak"),c(2)])
	
	model<-lm(get(ct) ~ Braak + Age + Sex , data = subset(predCCBest,Dataset == "LBB1"))
	resOut<-rbind(resOut, c("Samples" = "LBB1", "CellType" = ct, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))		
	te<-c(te, summary(model)$coefficients[c("Braak"),c(1)])
	seTE<-c(seTE, summary(model)$coefficients[c("Braak"),c(2)])
	
	model<-lm(get(ct) ~ Braak + Age + Sex , data = subset(predCCBest,Dataset == "MS"))
	resOut<-rbind(resOut, c("Samples" = "MS", "CellType" = ct, t(summary(model)$coefficients[c("Braak", "Age", "Sex2"),c(1,2,4)])))		
	te<-c(te, summary(model)$coefficients[c("Braak"),c(1)])
	seTE<-c(seTE, summary(model)$coefficients[c("Braak"),c(2)])
	
	## meta analysis
	outMeta<-metagen(te, seTE)
	resMeta<-rbind(resMeta, c(ct,outMeta$TE.fixed,outMeta$seTE.fixed,outMeta$pval.fixed, outMeta$TE.random, outMeta$seTE.random,outMeta$pval.random, outMeta$tau, outMeta$I2, outMeta$Q,1-pchisq(outMeta$Q, outMeta$df.Q)))
	
	reg<-lm(get(ct) ~ Braak, data = predCCBest)
	
	fig1a[[j]]<-ggplot(subset(predCCBest, !is.na(Braak)), aes_string(x="Braak_factor", y=ct, fill = "Braak_factor"))  +
	  geom_violin(scale = 'width')  +
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos, col = "white") + theme(legend.position = "none") + ylim(0,1) + geom_abline(intercept = coef(reg)[1], slope = coef(reg)[2])
				   
	
}	

ggarrange(plotlist=fig1a, nrow = 2, ncol = 5)
ggsave(filename = file.path(plotPath, "ViolinPlotCellCompBestByBraak.pdf"),  units = "in", width = 18, height = 8)


write.csv(resOut, file.path(plotPath, "CellCompAgainstBraakStatusByDataset.csv"))
write.csv(resMeta, file.path(plotPath, "CellCompAgainstBraakStatusMetaAnalysis.csv"))

## define AD status
predCCBest$AD<-rep(NA, nrow(predCCBest))
predCCBest$AD[which(predCCBest[,"Braak"] <3)]<-0
predCCBest$AD[which(predCCBest[,"Braak"] >4)]<-1
predCCBest$AD <- as.factor(predCCBest$AD)

resOut<-NULL
resMeta<-NULL
fig1a<-list()

for(j in 1:10){
	ct<-colnames(predCCBest)[j]

	te<-NULL
	seTE<-NULL

	model<-lm(get(ct) ~ AD + Age + Sex + Dataset + AD*Dataset, data = predCCBest)
	resOut<-rbind(resOut, c("Samples" = "All","CellType" = ct, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))
	
	model<-lm(get(ct) ~ AD + Age + Sex , data = subset(predCCBest,Dataset == "BDR"))
	resOut<-rbind(resOut, c("Samples" = "BDR", "CellType" = ct, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))		
	te<-c(te, summary(model)$coefficients[c("AD1"),c(1)])
	seTE<-c(seTE, summary(model)$coefficients[c("AD1"),c(2)])
	
	model<-lm(get(ct) ~ AD + Age + Sex , data = subset(predCCBest,Dataset == "LBB1"))
	resOut<-rbind(resOut, c("Samples" = "LBB1", "CellType" = ct, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))		
	te<-c(te, summary(model)$coefficients[c("AD1"),c(1)])
	seTE<-c(seTE, summary(model)$coefficients[c("AD1"),c(2)])
	
	model<-lm(get(ct) ~ AD + Age + Sex , data = subset(predCCBest,Dataset == "MS"))
	resOut<-rbind(resOut, c("Samples" = "MS", "CellType" = ct, t(summary(model)$coefficients[c("AD1", "Age", "Sex2"),c(1,2,4)])))		
	te<-c(te, summary(model)$coefficients[c("AD1"),c(1)])
	seTE<-c(seTE, summary(model)$coefficients[c("AD1"),c(2)])
	
	## meta analysis
	outMeta<-metagen(te, seTE)
	resMeta<-rbind(resMeta, c(ct,outMeta$TE.fixed,outMeta$seTE.fixed,outMeta$pval.fixed, outMeta$TE.random, outMeta$seTE.random,outMeta$pval.random, outMeta$tau, outMeta$I2, outMeta$Q,1-pchisq(outMeta$Q, outMeta$df.Q)))
	
	
	fig1a[[j]]<-ggplot(subset(predCCBest, !is.na(AD)), aes_string(x="AD", y=ct, fill = "AD"))  +
	  geom_violin(scale = 'width')  +
	  stat_summary(fun = "mean", 
				   geom = "point", 
				   position = pos, col = "white") + theme(legend.position = "none") + ylim(0,1) 
				   
	
}	

ggarrange(plotlist=fig1a, nrow = 2, ncol = 5)
ggsave(filename = file.path(plotPath, "ViolinPlotCellCompBestByAD.pdf"),  units = "in", width = 18, height = 8)


write.csv(resOut, file.path(plotPath, "CellCompAgainstADStatusByDataset.csv"))
write.csv(resMeta, file.path(plotPath, "CellCompAgainstADStatusMetaAnalysis.csv"))
