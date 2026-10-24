##---------------------------------------------------------------------#
##
## Title: Annotate EPIC manifest with gene lists
##
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# script adapted with thanks from:
# https://github.com/alicemfr/DevCortexDNAm/blob/main/enrichmentAnalyses/Genelist_enrichment/annotEPICman_genelists.r

#----------------------------------------------------------------------#
# Run 
#----------------------------------------------------------------------#

## Annotate EPIC manifest with gene lists

library(data.table)
library(tidyr)

refPath <- "references/EPICArray"

uniqueAnno <- function(row){ if(is.na(row)){row=''}; if(row != ""){ return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";")) } else { return(row) } }

setwd("Schizophrenia/4_analysis/")


#1. Load gene list files ========================================================================================================

# sczGWAS genes ------------------------------------------------------------------------------------------------------------
sczGWAS <- read.csv("SZ_GWAS_2022_prioritised.csv") # 685
sczGWAS <- sczGWAS[which(sczGWAS$Prioritised==1),] # limit to sczGWAS genes where priority == 1. n=120
#table(sczGWAS$syndromic) # 0=118; 1=115
sczGWAS <- sczGWAS[,c("Index.SNP", "Symbol.ID" ,"Ensembl.ID"  )]
colnames(sczGWAS)[1:2] <- c('SNP','Gene')

write.csv(sczGWAS, file="sczGWAS_genes_prioritised.csv")

##### genes which aren't working... 
# 1 - AC068490.2 
#  17 "CTD-2008L17.2"


# Schizophrenia SCHEMA genes ----------------------------------------------------------------------------------------------------
#scz <- read.csv("Singh2022_SCZ_GWAS.csv") # Supplementary Table 5 of Singh et al. (2022) https://doi.org/10.1038/s41586-022-04556-w
#schema <- scz[which(scz$Q.meta<0.05),]
#schema <- schema[,c('Gene.Symbol','Q.meta')]
#colnames(schema) <- c('Gene','Q-value')
#write.csv(schema, file="SCHEMA_genes.csv")

schema <- read.csv("SCHEMA_genes.csv", row.names=1)

# Combine -----------------------------------------------------------------------------------------------------------------------
# All genes across sczGWAS and SCHEMA. Duplicate genes present
df <- data.frame(Gene=c(sczGWAS$Gene, schema$Gene), GeneList=c(rep('sczGWAS',nrow(sczGWAS)),rep('SCHEMA',nrow(schema))))
write.csv(df, file="GeneList_sczGWAS_SCHEMA.csv")


#2. Load EPIC manifest ==========================================================================================================
epicManifest <- read.table(file.path(refPath, "EPIC.anno.GRCh38.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
epic <- epicManifest[,c("probeID","CHR","start","GeneNames","GeneClasses")] # limit to relevant columns
colnames(epic) <- c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")
epic$Gene <- unlist(lapply(epic$UCSC_RefGene_Name, uniqueAnno)) # extract unique transcripts of genes


#3. Annotate EPIC manifest with gene lists ======================================================================================

glist_name <- 'sczGWAS'
genelist <- sczGWAS
annotEPIC_glist <- function(epic, genelist, glist_name){
  epic[,glist_name] <- rep(0,nrow(epic))
  epic[,paste(glist_name,'Gene',sep='.')] <- rep('',nrow(epic))
  
  prg <- txtProgressBar()
  for(i in 2:nrow(genelist)){
    setTxtProgressBar(prg, i/nrow(genelist))

    print(i)
    
    # first try original gene name
    g <- as.character(genelist$Gene[i])
    print(g)
    indx <- grep( paste( paste0("^",g,"$"), paste0("^",g,";"), paste0(";",g,"$"), paste0(";",g,";"), sep='|'), epic$Gene) #regex = ^gene$|^gene;|;gene$|;gene;
    
    if(length(indx)>0){
      epic[indx,paste(glist_name,'Gene',sep='.')] <- g
    }
    
    # if not found, try alternative gene name
    if(length(indx)==0){ 
      print(paste(g, "had no matches."))
    }
      #a <- as.character(genelist$altGeneName[i])
      #if(g!=a){
       # print(paste("Trying alternative gene name:", a))
        #indx <- grep( paste( paste0("^",a,"$"), paste0("^",a,";"), paste0(";",a,"$"), paste0(";",a,";"), sep='|'), epic$Gene)
        #if(length(indx)!=0){
         # epic[indx,paste(glist_name,'Gene',sep='.')] <- a			
        #}else{
         # print(paste(a, "had no matches."))
        #}
      #}
    #}
    
    # set genelist column to 1 for probes annotated to genelist genes
    epic[indx,glist_name] <- 1
    
  }
  return(epic)
}

epic.annot <- epic

# Annotate EPIC manifest with sczGWAS genes...
epic.annot <- annotEPIC_glist(epic.annot, genelist=sczGWAS, glist_name='sczGWAS')

# ...and then with SCHEMA genes
epic.annot <- annotEPIC_glist(epic.annot, genelist=schema, glist_name='SCHEMA')



write.csv(epic.annot, file="EPIC_annot_sczGWAS_SCHEMAHG38.csv")
