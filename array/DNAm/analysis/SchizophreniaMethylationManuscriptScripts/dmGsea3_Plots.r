##---------------------------------------------------------------------#
##
## Title:dmGsea results plotting
##
##
## Input files: filtered results files from running dmGseaAnalysis.R script
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#


library(dplyr)
library(rrvgo)
library(GO.db)
library(igraph)
library(ggraph)
library(ggplot2)
library(org.Hs.eg.db)


dataDir <- "Schizophrenia"
outDir <- file.path(dataDir, "4_analysis/dmGseaV3") # main results
#outDir <- file.path(dataDir, "4_analysis/dmGseaHypoHyperAnalysis")
#outDir <- file.path(dataDir, "4_analysis/dmGseaNewThresholds")


setwd(outDir)


#----------------------------------------------------------------------#
# PLOT FUNCTION
#----------------------------------------------------------------------#

# https://github.com/ssayols/rrvgo/blob/master/vignettes/rrvgo.Rmd


makePlots <- function(inputfile, runType, GOtype, topN=50, thresh=500){
  
  setwd(file.path(outDir, runType))
  
  
  df <- sigRes[[inputfile]]
  df$GO_ID <- sub("\\|.*$", "", df$ID)
  
  df <- df[df$nGene < thresh, ] 
  
  simMatrix <- calculateSimMatrix(df$GO_ID,
                                  orgdb="org.Hs.eg.db",
                                  ont=GOtype,
                                  method="Rel")
  
  
  
  
  
  # replace zero FDRs before calculating scores if using gsProbe
  if(runType == "gsProbe"){
    df$FDR[df$FDR == 0] <- min(df$FDR[df$FDR > 0]) / 2
    scores <- setNames(-log10(df$FDR), df$GO_ID)
  } else {
    scores <- setNames(-log10(df$FDR), df$GO_ID)
  }
  
  
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  
  pdf(paste0("filteredResults/plots/", inputfile, "_RevigoPlots.pdf"))
  scatterPlot(simMatrix, reducedTerms)
  treemapPlot(reducedTerms)
  wordcloudPlot(reducedTerms, min.freq=1, colors="black")
  dev.off()
  
  
  
  # revigo top 50
  
  topReduced <- reducedTerms %>%
    arrange(desc(score))  %>%
    slice_head(n = topN)
  
  topReduced <- topReduced %>% dplyr::filter(go %in% rownames(simMatrix))
  
  pdf(paste0("filteredResults/plots/", inputfile, "_RevigoPlots_top", topN, ".pdf"))
  
  scatterPlot(simMatrix, topReduced)
  treemapPlot(topReduced)
  wordcloudPlot(topReduced, min.freq=1, colors="black")
  dev.off()
  
  
  #----------------------------------------------------------------------#
  # bubble Plots
  #----------------------------------------------------------------------#
  
  # all
  
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
  
  
  # top50
  
  topDF <- reducedTerms %>% arrange(desc(score)) %>% slice_head(n = topN)
  
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
  
  pdf(paste0("filteredResults/plots/", inputfile, "_plots_Top50BubblePlot.pdf"), height = 10, width = 20)
  print(p)
  dev.off()
  
  
  
  #----------------------------------------------------------------------#
  # network Plots
  #----------------------------------------------------------------------#
  
  # all
  
  topDF <- reducedTerms %>% arrange(desc(score))
  topGO <- topDF$go
  
  # Step 2: subset similarity matrix for top GO terms
  simMatrixTop <- simMatrix[topGO, topGO]
  
  # Step 3: create edge list for edges above a similarity threshold
  threshold <- 0.5  # only show edges where similarity > 0.5
  edges <- which(simMatrixTop > threshold & upper.tri(simMatrixTop), arr.ind = TRUE)
  edgeList <- data.frame(
    from = rownames(simMatrixTop)[edges[,1]],
    to   = colnames(simMatrixTop)[edges[,2]],
    weight = simMatrixTop[edges]
  )
  
  # Step 4: create
  g <- graph_from_data_frame(edgeList, vertices = topDF, directed = FALSE)
  
  # Step 5: plot network
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
  
  
  # Step 1: subset top 50
  
  topDF <- reducedTerms %>% arrange(desc(score)) %>% slice_head(n = topN)
  topGO <- topDF$go
  
  # Step 2: subset similarity matrix for top GO terms
  simMatrixTop <- simMatrix[topGO, topGO]
  
  # Step 3: create edge list for edges above a similarity threshold
  threshold <- 0.5  # only show edges where similarity > 0.5
  edges <- which(simMatrixTop > threshold & upper.tri(simMatrixTop), arr.ind = TRUE)
  edgeList <- data.frame(
    from = rownames(simMatrixTop)[edges[,1]],
    to   = colnames(simMatrixTop)[edges[,2]],
    weight = simMatrixTop[edges]
  )
  
  # Step 4: create
  g <- graph_from_data_frame(edgeList, vertices = topDF, directed = FALSE)
  
  # Step 5: plot network
  set.seed(42)
  netplot <- ggraph(g, layout = "fr") +   # Fruchterman-Reingold layout
    geom_edge_link(aes(width = weight), color = "grey70") +
    geom_node_point(aes(size = size, color = score)) +
    geom_node_text(aes(label = term), repel = TRUE, size = 3) +
    scale_color_viridis_c(option = "C") +
    theme_void() +
    labs(size = "Gene count", color = "-log10(FDR)", title = "GO Term Network")
  
  pdf(paste0("filteredResults/plots/", inputfile, "networkPlot_top50.pdf"))
  print(netplot)
  dev.off()
  
  
}




#----------------------------------------------------------------------#
# USE FUNCTION
#----------------------------------------------------------------------#

setwd(file.path(outDir, "gsGene"))


# fileList <- "filteredResults/gsGeneThreshold_GO_BP_nGeneMoreThan1.csv" # main results
fileList <- list.files("filteredResults/", "BP_nGeneMoreThan1.csv", full.names = T) # GO results only
sigRes<- lapply(fileList, function(x) read.csv(x, stringsAsFactors = F))
names(sigRes) <- tools::file_path_sans_ext(basename(fileList))


# loop over all files to make plots
for(name in names(sigRes)){
  
  print(name)
  ont <- sub(".*GO_([A-Z]{2}).*", "\\1", name) # extract BP/CC/MF
  print(ont)
  makePlots(name, "gsGene", ont)
  
}

