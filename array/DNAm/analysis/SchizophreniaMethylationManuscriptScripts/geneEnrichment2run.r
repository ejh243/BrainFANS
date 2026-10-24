##---------------------------------------------------------------------#
##
## Title: Gene list enrichments
##
## Purpose of script: perform gene list enrichment analysis of genes 
##						from the SCZ GWAS and SCHEMA gene lists
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# script adapted with thanks from:
# https://github.com/alicemfr/DevCortexDNAm/blob/main/enrichmentAnalyses/Genelist_enrichment/EWAS_genelistEnrichment.r


#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#

dataDir <- "Schizophrenia"
setwd(dataDir)


# load file created with script from Alice - https://github.com/alicemfr/DevCortexDNAm/blob/main/enrichmentAnalyses/Genelist_enrichment/annotEPICman_genelists.r
epicAnnotGeneList <- read.csv(file.path(dataDir, "4_analysis/EPIC_annot_sczGWAS_SCHEMAHG38.csv"), row.names=1)
sczGWAS <- unique(epicAnnotGeneList$sczGWAS.Gene[-which(is.na(epicAnnotGeneList$sczGWAS.Gene)|epicAnnotGeneList$sczGWAS.Gene=='')])
schema <- unique(epicAnnotGeneList$SCHEMA.Gene[-which(is.na(epicAnnotGeneList$SCHEMA.Gene)|epicAnnotGeneList$SCHEMA.Gene=='')])
glistFile <- data.frame(Gene=c(sczGWAS, schema), GeneList=c(rep('sczGWAS',length(sczGWAS)),rep('SCHEMA',length(schema))))
glistFile <- rbind(glistFile, data.frame(Gene=unique(glistFile$Gene), GeneList=rep('Combined', length(unique(glistFile$Gene)))))


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTIONS
#----------------------------------------------------------------------#

apply_tests<-function(glist, pathway, gene_test, gene_size){

	model <- glm(pathway ~ gene_test + gene_size, family="binomial")
	p1 <- summary(model)$coefficients['gene_test','Pr(>|z|)']
	p2 <- summary(model)$coefficients['gene_size','Pr(>|z|)']
	b1 <- summary(model)$coefficients[c("gene_test"),'Estimate']
	b2 <- summary(model)$coefficients[c("gene_size"),'Estimate']
	or <- exp(b1)
	
	# divide the p-value by 2 to convert two-sided to one-sided test
	if(or > 1){															                        
		p1 <- p1/2
	}else{
		p1 <- 1-(p1/2)
	}

		  # Find gene names in pathway that have >= 1 DMP
  gene_names <- names(gene_test)[which(pathway == 1 & gene_test >= 1)]
  gene_names_collapsed <- paste(gene_names, collapse = ";")

	
	return(c(
	glist,                                                    # the genelist we've tested
	sum(gene_size[which(pathway == 1)]),                      # total number of probes annotated to genelist genes
	length(which(pathway == 1)),                              # number of genes in genelist
	sum(gene_test[which(pathway == 1  & gene_test >= 1)]),    # total number of DMPs annotated to genelist genes
	length(which(pathway == 1 & gene_test >= 1)),             # number of genelist genes that contain a DMP
	b1,                                                       # effect size for gene_test
	or,                                                       # odds ratio for gene_test
	p1,                                                       # p-value for gene_test
	b2,                                                       # effect size for gene_size
	p2,
	gene_names_collapsed                                        # gene names in pathway with >=1 DMP                                                       # p-value for gene_size
	))
	

}

geneEnrichment <- function(inputFile, resultCol, epicAnnotGeneList, outputFile, glists, glistFile, pThresh){

	print(paste('Input file:', inputFile))
	
	if (grepl(".rdat", inputFile)) {
    	loaded_obj <- load(inputFile)  # This returns the name(s) of the object(s) loaded
   	 res <- as.data.frame(get(loaded_obj))  # Access the loaded object dynamically
	}
	
	print(dim(res))
    res$IlmnID <- row.names(res)

	# annotate EWAS results with EPIC manifest that has been pre-annotated with whether probe is in genelist
	epicMan <- epicAnnotGeneList[match(rownames(res), epicAnnotGeneList$IlmnID),]
	res <- cbind(res, as.data.frame(epicMan))
	if(any(res$CHR %in% 'Y')){ res <- res[-which(res$CHR=='Y'),] }
	print(dim(res))


	#1. Number of probes per gene
	gene_size <- table(unlist(strsplit(res$Gene, "\\;")))
	print(length(gene_size))
	
	#2. Number of DMPs annotated to gene
	res$DMP <- as.numeric(res[,resultCol]<pThresh)
	gene_dmp <- table(unlist(strsplit(res$Gene[res$DMP==1], "\\;")))
	gene_test <- gene_dmp[match(names(gene_size), names(gene_dmp))]
	names(gene_test) <- names(gene_size)
	gene_test[is.na(gene_test)]<-0
	print(length(gene_test))
	
	mat <- matrix(NA,nrow=length(glists),ncol=11)
	colnames(mat) <- c("GeneList", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "Estimate.GenesinTestList", "OR.GenesinTestList", "P.GenesinTestList", "Estimate.GenesSize", "P.GeneSize", "GeneNamesInPathway")
	
	for(i in 1:length(glists)){
	  glist <- glists[i]
	  print(glist)
	  genes <- glistFile[glistFile$GeneList==glist,'Gene']
	  
	  #3. Is gene in genelist: 0/1
	  pathway <- as.numeric(names(gene_size) %in% genes)
	  print(length(pathway))
	  
	  mat[i,] <- apply_tests(glist, pathway, gene_test, gene_size)
	  
	}
	
	write.csv(mat, outputFile)


# extract dmps
for (glist in glists) {
  genes <- glistFile[glistFile$GeneList == glist, 'Gene']
  
  # DMPs in this gene list
  dmps_in_list <- res[res$DMP == 1 & sapply(strsplit(res$Gene, ";"), function(x) any(x %in% genes)), ]
  
  # Save a separate CSV per gene list
  dmpFile <- paste0(tools::file_path_sans_ext(outputFile), "_DMPs_", glist, ".csv")
  write.csv(dmps_in_list, dmpFile, row.names = FALSE)
  
  print(paste("DMP table for", glist, "saved to:", dmpFile))
  }

}


#----------------------------------------------------------------------#
# RUN ANALYSIS
#----------------------------------------------------------------------#

geneEnrichment(inputFile = "4_analysis/EWAS/new/NeuN+LM_FDR.rdata", resultCol = "NullModel_SCZ_P", epicAnnotGeneList = epicAnnotGeneList, outputFile = "geneEnrichment_NeuN+LM_1e04HG38.csv", glists=c('sczGWAS','SCHEMA', 'Combined'), glistFile = glistFile, pThresh = 1e-4)

