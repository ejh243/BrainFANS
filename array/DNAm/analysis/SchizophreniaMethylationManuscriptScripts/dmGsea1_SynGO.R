##---------------------------------------------------------------------#
##
## Title:dmGsea
##
##
## Input files: EWAS results file (neun only)
##              betas matrix
##              
##
## Output files: results tables for gsGene and gsProbe
##---------------------------------------------------------------------#


# User defined geneset


#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# https://bioconductor.org/packages/devel/bioc/vignettes/dmGsea/inst/doc/dmGsea.html


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(dmGsea)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(dplyr)

# plotting
library(rrvgo)
library(GO.db)
library(igraph)
library(ggraph)


dataDir <- "Schizophrenia"
normData <- file.path(dataDir, "3_normalised/normalised.rdata")
resData <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")


outDir <- file.path(dataDir, "4_analysis/dmGseaV3/syngo1.3_complete_data/results")

source(file.path(dataDir, "scripts/dmGsea_gsGene.R"))
source(file.path(dataDir, "scripts/dmGsea_utils.R"))
source(file.path(dataDir, "scripts/dmGsea_gsRank.R"))
source(file.path(dataDir, "scripts/dmGsea_gsProbe.R"))
source(file.path(dataDir, "scripts/dmGsea_gsPG.R"))



#----------------------------------------------------------------------#
# LOAD AND FORMAT DATA
#----------------------------------------------------------------------#

load(resData)
neun <- res[[2]]
neun$Name <- row.names(neun)

load(normData)

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to cell types used
cellTypes <- "NeuN+"
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% cellTypes,]

# subset betas matrix to same samples
celltypeNormbeta <- celltypeNormbeta[, QCmetrics$Basename]
celltypeNormbeta <- celltypeNormbeta[complete.cases(celltypeNormbeta),]
celltypeNormbeta <- celltypeNormbeta[rownames(neun),]


# create probe.p object
neun.p <- neun %>% 
  dplyr::select(Name, NullModel_SCZ_P) %>%
  dplyr::rename(p = NullModel_SCZ_P)



#----------------------------------------------------------------------#
# Create GeneSet
#----------------------------------------------------------------------#

# load synapse GO data file, and original synGO data file for entrezIDs
nickDat <- read.csv(file.path(dataDir, "4_analysis/dmGseaV3/syngo1.3_complete_data/ontologies.csv"))
nickTab <- read.table(file.path(dataDir, "4_analysis/dmGseaV3/syngo1.3_complete_data/SynGO-2023_CC-BP_ensembl_FromNick.txt"))


# reformat name so can match with nickTab
nickDat$nameReformat <- gsub(" ", "_", nickDat$name)
nickDat$nameReformat <- gsub("\\(", "", nickDat$nameReformat)
nickDat$nameReformat <- gsub("\\)", "", nickDat$nameReformat)
nickDat$nameReformat <- gsub(":", "", nickDat$nameReformat)
nickDat$nameReformat <- gsub(",", "", nickDat$nameReformat)


# format to match example above
toUse <- nickDat[nickDat$nameReformat %in% nickTab$V1,]

toUseBP <- toUse[toUse$domain == "BP",]
toUseCC <- toUse[toUse$domain == "CC",]


# all together 
userGeneset <- tapply(toUse$entrez_id, toUse$nameReformat, function(x) {
  trimws(unlist(strsplit(as.character(x), ",")))
})


# BP and CC seperately 
userGenesetBP <- tapply(toUseBP$entrez_id, toUseBP$nameReformat, function(x) {
  trimws(unlist(strsplit(as.character(x), ",")))
})

userGenesetCC <- tapply(toUseCC$entrez_id, toUseCC$nameReformat, function(x) {
  trimws(unlist(strsplit(as.character(x), ",")))
})



#----------------------------------------------------------------------#
# Run
#----------------------------------------------------------------------#

setwd(outDir)


# BP and CC together
gsGeneLocal(probe.p = neun.p, 
            Data4Cor=celltypeNormbeta, 
            arrayType="EPIC", 
            ncore=1,
            method="Threshold",
            geneSet=userGeneset,
            outfile="gsGeneThresholdSynapseGO",
            outGenep=TRUE)


# BP only
gsGeneLocal(probe.p = neun.p, 
            Data4Cor=celltypeNormbeta, 
            arrayType="EPIC", 
            ncore=1,
            method="Threshold",
            geneSet=userGenesetBP,
            outfile="gsGeneThresholdSynapseGO_BP",
            outGenep=FALSE)


# CC only
gsGeneLocal(probe.p = neun.p, 
            Data4Cor=celltypeNormbeta, 
            arrayType="EPIC", 
            ncore=1,
            method="Threshold",
            geneSet=userGenesetCC,
            outfile="gsGeneThresholdSynapseGO_CC",
            outGenep=FALSE)






#----------------------------------------------------------------------#
# Filter to just sig results
#----------------------------------------------------------------------#

gsRes <- read.csv("gsGeneThresholdSynapseGO_userSet_userSet.csv")


sigRes <- gsRes %>% filter(nGene > 1 & FDR < 0.05) %>%
    arrange(FDR) 

write.csv(sigRes, "gsGeneThresholdSynapseGO_Filtered.csv")



# read in results files
fileList <- list.files(".", pattern = "gsGeneThresholdSynapseGO", full.names = T)
resFiles<- lapply(fileList, function(x) read.csv(x, stringsAsFactors = F))
names(resFiles) <- tools::file_path_sans_ext(basename(fileList))


# subset to jsut sig results
sigRes <- lapply(resFiles, function(df) {
  df %>%   filter(nGene > 1 & FDR < 0.05) %>%
    arrange(FDR) %>%
    mutate(GO_ID = sub(".*_(GO[0-9]+)$", "\\1", ID))
})


# save filtered results
for(name in names(sigRes)){
  outFile <- paste0("filteredResults/", name, "_nGeneMoreThan1.csv")
  write.csv(sigRes[[name]], outFile)
}


#----------------------------------------------------------------------#
# plot -  only for BP
#----------------------------------------------------------------------#

# not filtered to rm pathwayw with >500 genes 
# removed anything doesn't have a GO ID

inputfile <- "gsGeneThresholdSynapseGO_BP_userSet_userSet"

df <- sigRes[[inputfile]]

# remove anything doesn't have a GO ID
df <- df[grepl("^GO", df$GO_ID), ]

# add colon after "GO
df$GO_ID <- sub("(GO)([0-9]+)", "\\1:\\2", df$GO_ID)


simMatrix <- calculateSimMatrix(df$GO_ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(df$FDR), df$GO_ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")


pdf(paste0("filteredResults/plots/", inputfile, "_RevigoPlots.pdf"))
scatterPlot(simMatrix, reducedTerms)
treemapPlot(reducedTerms)
wordcloudPlot(reducedTerms, min.freq=1, colors="black")
dev.off()


# bubble plot
topDF <- reducedTerms %>% arrange(desc(score))

p <- ggplot(topDF, aes(x = reorder(term, score), y = score, size = size, color = score)) +
  geom_point() +
  coord_flip() +
  scale_color_viridis_c() +
  labs(
    x = "GO Term",
    y = "-log10(FDR)",
    size = "Gene count",
    color = "Significance"
  ) +
  theme_minimal()

pdf(paste0("filteredResults/plots/", inputfile, "_plots_BubblePlot.pdf"), height = 10, width = 20)
print(p)
dev.off()



topDF <- reducedTerms %>% arrange(desc(score))
topGO <- topDF$go

# subset similarity matrix for top GO terms
simMatrixTop <- simMatrix[topGO, topGO]

# Step 3: create edge list for edges above a similarity threshold
threshold <- 0.5  # only show edges where similarity > 0.5
edges <- which(simMatrixTop > threshold & upper.tri(simMatrixTop), arr.ind = TRUE)
edgeList <- data.frame(
  from = rownames(simMatrixTop)[edges[,1]],
  to   = colnames(simMatrixTop)[edges[,2]],
  weight = simMatrixTop[edges]
)

# create
g <- graph_from_data_frame(edgeList, vertices = topDF, directed = FALSE)

# plot network
set.seed(42)
netplot <- ggraph(g, layout = "fr") +   # Fruchterman-Reingold layout
  geom_edge_link(aes(width = weight), color = "grey70") +
  geom_node_point(aes(size = size, color = score)) +
  geom_node_text(aes(label = term), repel = TRUE, size = 3) +
  scale_color_viridis_c(option = "C") +
  theme_void() +
  labs(size = "Gene count", color = "-log10(FDR)", title = "GO Term Network")

pdf(paste0("filteredResults/plots/", inputfile, "networkPlot.pdf"))
print(netplot)
dev.off()



