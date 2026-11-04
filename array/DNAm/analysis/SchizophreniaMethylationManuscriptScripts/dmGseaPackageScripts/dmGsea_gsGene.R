##
## Title: dmGsea gsGene.R functions
##
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# Running the vignette code from here gives error: "Diagonal values in 'R' must all be equal to 1."
# https://bioconductor.org/packages/release/bioc/vignettes/dmGsea/inst/doc/dmGsea.html#example-1-differentially-methylated-probes-from-ewas

# Looking at the code it appears to be caused when there is a discrepancy between CpG IDs 
# in the betas matrix/results and the manifest file (e.g. probes removed during QC)
# to fix this the following has been added / modified in the combinep function:
#
#  common <- subset$Name[subset$Name %in% rownames(Data4Cor)]
#  dat <- t(Data4Cor[common, , drop=FALSE])
#  #dat <- t(Data4Cor[as.vector(subset$Name),]) - this one replaced with the two above.
#
# This ensures that that the dat object is only created using CpGs that in the dataset.
#
# In order to run the main gsGene function, the following functions also need to be overwritten 
# as the package is 'locked' to overwriting the combinep function.



#----------------------------------------------------------------------#
# EDITED FUNCTIONS
#----------------------------------------------------------------------#

# this is the function with edited code
combinepLocal <- function(gid,Data4Cor,GeneProbeTablep,combpMethod,combpAdjust){
  #    if(!requireNamespace("poolr",quietly=TRUE))stop("poolr package required")
  subset <- GeneProbeTablep[GeneProbeTablep$entrezid==gid,]
  allp <- subset$p
  if(length(allp)>=2){
    if(is.null(Data4Cor)){
      eval(parse(text=paste0("poolr::",combpMethod,"(p=allp,adjust='none')")))$p
    }else{
      common <- subset$Name[subset$Name %in% rownames(Data4Cor)]
      dat <- t(Data4Cor[common, , drop=FALSE])
      # dat <- t(Data4Cor[as.vector(subset$Name),])
      Rmatrix <- cor(dat,use="pairwise.complete.obs")
      diag(Rmatrix) <- rep(1, nrow(Rmatrix))
      eval(parse(text=paste0("poolr::",combpMethod,"(p=allp,adjust=combpAdjust,
            R=Rmatrix)")))$p
    }}else{allp}
}





# edited to run combinepLocal
calcGenePLocal <- function(Data4Cor, GeneProbeTable, combpMethod,
                           combpAdjust, ncore, outGenep, species) {
  genelist <- unique(GeneProbeTable$entrezid)
  genep <- mclapply(genelist, combinepLocal, Data4Cor = Data4Cor,
                    GeneProbeTablep = GeneProbeTable, combpMethod = combpMethod,
                    combpAdjust = combpAdjust, mc.cores = ncore)
  genep <- unlist(genep)
  names(genep) <- genelist
  genep <- genep[!is.na(genep)]
  if (outGenep) {
    tmpgn <- geneID2geneName(names(genep), species)
    tmp <- data.frame(geneid = names(genep), gene = tmpgn, p = genep)
    write.csv(tmp, "gsGene_genep.csv", row.names = FALSE)
  }
  genep
}




# edited to run calcGenePLocal
gsGeneLocal <- function(probe.p, Data4Cor = NULL, method = "Threshold",
                        FDRthre = 0.05, nTopGene = NULL, GeneProbeTable = NULL,
                        arrayType = NULL, gSetName = "KEGG", geneSet = NULL,
                        species = "Human", combpMethod = "fisher",
                        combpAdjust = "nyholt", outfile = "gsGene", outGenep = FALSE,
                        gseaParam = 1, nperm = 1e4, ncore=1) {
  if (method == "Ranking") nTopGene <- NULL
  if (is.null(Data4Cor)) combpAdjust <- "none"
  input <- processInput(probe.p, Data4Cor, method, arrayType, GeneProbeTable,
                        combpAdjust)
  Data4Cor <- input$Data4Cor
  GeneProbeTable <- input$GeneProbeTable
  probe.p <- input$probe.p
  genelist <- unique(GeneProbeTable$entrezid)
  genep <- calcGenePLocal(Data4Cor, GeneProbeTable, combpMethod,
                          combpAdjust, ncore, outGenep, species)
  gSetName12 <- ParseGsetName(geneSet, gSetName)
  for (gset1 in gSetName12$gSetName1) {
    for (gset2 in gSetName12$gSetName2) {
      geneSet <- processGeneSets(geneSet, gset1, gset2, species,
                                 unique(names(genep)))
      if (method == "Ranking") {
        resu <- processRankingMethod(geneSet, genep, nperm, ncore, gseaParam)
      } else if (method == "Threshold") {
        resu <- processThresholdMethod(geneSet,genep,nTopGene,FDRthre,
                                       species)
      } else {
        stop("method need to be Threshold or Ranking")
      }
      write.csv(resu, paste0(outfile, "_", gset1, "_", gset2, ".csv"),
                row.names = FALSE)
    }
  }
}




#----------------------------------------------------------------------#
# C+P Functions
#----------------------------------------------------------------------#


# permutation test simulation
permuTestSim <- function(rgl, nSim, ncore, maxGsSize, gsSize) {
  ng <- length(rgl)
  resu <- mclapply(seq_len(nSim), function(k) {
    rgs <- sample.int(ng, maxGsSize)
    vapply(gsSize, function(gsSize_i) {
      rgs.i <- sort(rgs[seq_len(gsSize_i)])
      rglsumi <- sum(rgl[rgs.i])
      ngi <- ng - gsSize_i
      x <- rgl[rgs.i]
      y <- rgs.i - c(seq_len(gsSize_i))
      max(cumsum(x) / rglsumi - y / ngi)
    }, FUN.VALUE = numeric(1))
  }, mc.cores = ncore)
  mat.null <- do.call(rbind, lapply(resu, unlist))
  colnames(mat.null) <- gsSize
  mat.null
}

# permutation P-value calculation
permuTestP <- function(mat.null, gsSize0, es0, nSim) {
  gsSize0 <- as.character(gsSize0)
  permuP <- vapply(seq_len(length(es0)), function(i) {
    mean(mat.null[, gsSize0[i]] >= es0[i])
  }, FUN.VALUE = numeric(1))
  permuP[permuP == 0] <- 1 / (nSim * 2)
  permuP
}

# beta fit P-value calculation
permuTestBeta <- function(mat.null, gsSize0, es0) {
  b_para <- do.call(rbind,lapply(seq_len(ncol(mat.null)), fitBeta, 
                                 dat = mat.null))
  rownames(b_para) <- colnames(mat.null)
  b_para <- b_para[as.character(as.vector(gsSize0)), ]
  p <- vapply(seq_len(length(es0)), function(i) {
    exp(pbeta(es0[i], shape1 = b_para[i, 1], shape2 = b_para[i, 2],
              log.p = TRUE, lower.tail = FALSE))
  }, FUN.VALUE = numeric(1))
  p
}

# permutation test function
permuTest2 <- function(gSet, rgl, nSim = 10000, es0, ncore) {
  gsSize0 <- vapply(gSet, length, FUN.VALUE = numeric(1))
  gsSize <- unique(gsSize0)
  maxGsSize <- max(gsSize)
  mat.null <- permuTestSim(rgl, nSim, ncore, maxGsSize, gsSize)
  permuP <- permuTestP(mat.null, gsSize0, es0, nSim)
  p <- permuTestBeta(mat.null, gsSize0, es0)
  return(list(p = p, permuP = permuP))
}

# Function to calculate gene set P-value
p.gs.gene <- function(gs, sig.gene, nGene, nSig) {
  gs.ng <- length(gs)
  gs.nSig <- sum(gs %in% sig.gene) - 1
  phyper(gs.nSig, gs.ng, nGene - gs.ng, nSig, lower.tail = FALSE)
}

# function to validate input data
processInput <- function(probe.p, Data4Cor, method, arrayType, GeneProbeTable,
                         combpAdjust){
  if (!is.null(Data4Cor)){
    if (is(Data4Cor, "SummarizedExperiment")) {
      Data4Cor <- assays(Data4Cor)$beta
    }
    if (!is(Data4Cor, "matrix")) {
      stop("Data4Cor should be a numeric matrix or SummarizedExperiment")
    }}
  if (is.null(GeneProbeTable)) {
    if (!(arrayType %in% c("450K", "EPIC"))) {
      stop("Specify arrayType as 450K or EPIC, or provide GeneProbeTable")
    }
    GeneProbeTable <- getIlluminaAnnotation(arrayType = arrayType)
  }
  
  names(probe.p)[which(names(probe.p) == "probe")] <- "Name"
  names(probe.p)[which(names(probe.p) == "P")] <- "p"
  probe.p$Name <- as.character(probe.p$Name)
  if (combpAdjust != "none") {
    tmp <- sum(!probe.p$Name %in% rownames(Data4Cor))
    if (tmp > 0) {
      probe.p <- probe.p[probe.p$Name %in% rownames(Data4Cor), ]
      message(tmp, " probes were missing in Data4Cor, excluded from analysis")
    }
  }
  
  GeneProbeTable$p <- probe.p$p[match(GeneProbeTable$Name, probe.p$Name)]
  GeneProbeTable <- GeneProbeTable[complete.cases(GeneProbeTable), ]
  GeneProbeTable <- GeneProbeTable[order(GeneProbeTable$entrezid,
                                         GeneProbeTable$p), ]
  GeneProbeTable$Name <- as.character(as.vector(GeneProbeTable$Name))
  GeneProbeTable$entrezid <- as.character(as.vector(GeneProbeTable$entrezid))
  GeneProbeTable$p <- as.numeric(as.vector(GeneProbeTable$p))
  list(Data4Cor = Data4Cor, GeneProbeTable = GeneProbeTable,probe.p=probe.p)
}

# calculate gene P-values
# calcGeneP <- function(Data4Cor, GeneProbeTable, combpMethod,
#                       combpAdjust, ncore, outGenep, species) {
#   genelist <- unique(GeneProbeTable$entrezid)
#   genep <- mclapply(genelist, combinep, Data4Cor = Data4Cor,
#                     GeneProbeTablep = GeneProbeTable, combpMethod = combpMethod,
#                     combpAdjust = combpAdjust, mc.cores = ncore)
#   genep <- unlist(genep)
#   names(genep) <- genelist
#   genep <- genep[!is.na(genep)]
#   if (outGenep) {
#     tmpgn <- geneID2geneName(names(genep), species)
#     tmp <- data.frame(geneid = names(genep), gene = tmpgn, p = genep)
#     write.csv(tmp, "gsGene_genep.csv", row.names = FALSE)
#   }
#   genep
# }

# process gene sets
processGeneSets <- function(geneSet, gset1, gset2, species, universe) {
  if (is.null(geneSet)) {
    if (gset1 == "KEGG") {
      geneSet <- getKEGG(species)
    } else if (gset1 == "GO") {
      geneSet <- getGO(gset2, species)
    } else if (gset1 == "MSigDB") {
      geneSet <- getMSigDB(gset2, species)
    } else if (gset1 == "Reactome") {
      geneSet <- getReactome(species)
    } else {
      stop("Specify gSetName as KEGG, GO, MSigDB, or Reactome, or provide
                geneSet")
    }
  }
  
  if (!is.list(geneSet)) stop("geneSet should be a list object")
  geneSet <- lapply(geneSet, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- x[x %in% universe]
    unique(x)
  })
  inUniv <- vapply(geneSet, function(x) length(x) > 0, FUN.VALUE = logical(1))
  geneSet[inUniv]
}

# ranking method
processRankingMethod <- function(geneSet, genep, nperm, ncore, gseaParam) {
  genep[genep == 0] <- min(genep[!genep == 0])
  genep[genep == 1] <- max(genep[!genep == 1])
  genep <- sort(-log(genep), decreasing = TRUE)
  genep <- genep^gseaParam
  es0 <- vapply(geneSet, function(gsi) max(es(gsi, genep)),
                FUN.VALUE = numeric(1))
  p <- permuTest2(gSet = geneSet, rgl = genep, nSim = nperm, es0 = es0,
                  ncore = ncore)
  nGene <- vapply(geneSet, length, FUN.VALUE = numeric(1))
  data.frame(ID = names(geneSet), nGene = nGene, EnrichScore = es0,
             p = p$p, permuP = p$permuP, FDR = p.adjust(p$p, method = "fdr"))
}

# Helper function for threshold method
processThresholdMethod <- function(geneSet, genep, nTopGene, FDRthre, species){
  all.gene <- unique(names(genep))
  if (!is.null(nTopGene)) {
    sig.gene <- names(genep)[order(genep)][seq_len(nTopGene)]
  } else {
    sig.gene <- names(genep)[p.adjust(genep, method = "fdr") < FDRthre]
  }
  if (length(sig.gene) < 10) {
    stop("Number of significant genes<10, adjust threshold or use ranking")
  }
  sgName <- geneID2geneName(sig.gene, species)
  nSigGene <- length(sig.gene)
  nGene <- length(all.gene)
  p <- vapply(geneSet, p.gs.gene, sig.gene = sig.gene, nGene = nGene,
              nSig = nSigGene, FUN.VALUE = numeric(1))
  nGene <- vapply(geneSet, length, FUN.VALUE = numeric(1))
  nSigGene <- vapply(geneSet, function(x) sum(sig.gene %in% x),
                     FUN.VALUE = numeric(1))
  Sig.Gene <- vapply(geneSet, 
                     function(x) paste(sgName[sig.gene[sig.gene %in% x]],collapse = ";"),
                     FUN.VALUE = character(1))
  data.frame(ID = names(geneSet), nGene = nGene, nSigGene = nSigGene, p = p,
             FDR = p.adjust(p, method = "fdr"), Sig.Gene = Sig.Gene)
}

# Main function
# gsGene <- function(probe.p, Data4Cor = NULL, method = "Threshold",
#                    FDRthre = 0.05, nTopGene = NULL, GeneProbeTable = NULL,
#                    arrayType = NULL, gSetName = "KEGG", geneSet = NULL,
#                    species = "Human", combpMethod = "fisher",
#                    combpAdjust = "nyholt", outfile = "gsGene", outGenep = FALSE,
#                    gseaParam = 1, nperm = 1e4, ncore=1) {
#   if (method == "Ranking") nTopGene <- NULL
#   if (is.null(Data4Cor)) combpAdjust <- "none"
#   input <- processInput(probe.p, Data4Cor, method, arrayType, GeneProbeTable,
#                         combpAdjust)
#   Data4Cor <- input$Data4Cor
#   GeneProbeTable <- input$GeneProbeTable
#   probe.p <- input$probe.p
#   genelist <- unique(GeneProbeTable$entrezid)
#   genep <- calcGeneP(Data4Cor, GeneProbeTable, combpMethod,
#                      combpAdjust, ncore, outGenep, species)
#   gSetName12 <- ParseGsetName(geneSet, gSetName)
#   for (gset1 in gSetName12$gSetName1) {
#     for (gset2 in gSetName12$gSetName2) {
#       geneSet <- processGeneSets(geneSet, gset1, gset2, species,
#                                  unique(names(genep)))
#       if (method == "Ranking") {
#         resu <- processRankingMethod(geneSet, genep, nperm, ncore, gseaParam)
#       } else if (method == "Threshold") {
#         resu <- processThresholdMethod(geneSet,genep,nTopGene,FDRthre,
#                                        species)
#       } else {
#         stop("method need to be Threshold or Ranking")
#       }
#       write.csv(resu, paste0(outfile, "_", gset1, "_", gset2, ".csv"),
#                 row.names = FALSE)
#     }
#   }
# }

