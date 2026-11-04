##
## Title: dmGsea utils.R functions
##
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# Running the vignette code from here gives error: "Diagonal values in 'R' must all be equal to 1."
# https://bioconductor.org/packages/release/bioc/vignettes/dmGsea/inst/doc/dmGsea.html#example-1-differentially-methylated-probes-from-ewas

# Looking at the code it appears to be caused when there is a discrepency between CpG IDs 
# in the betas matrix/results and the manifest file (e.g. probes removed during QC)
# to fix this the following has been added / modifided in the combinep function:
#
#  common <- subset$Name[subset$Name %in% rownames(Data4Cor)]
#  dat <- t(Data4Cor[common, , drop=FALSE])
#  #dat <- t(Data4Cor[as.vector(subset$Name),]) - this one replaced with the two above.
#
# This ensures that that the dat object is only created using CpGs that in the dataset.
#
# In order to run the main gsGene function, the following fucntions also need to be overwritten 
# as the package is 'locked' to overwritting the combinep function.



#----------------------------------------------------------------------#
# EDITED FUNCTIONS
#----------------------------------------------------------------------#

# None


#----------------------------------------------------------------------#
# C+P Functions
#----------------------------------------------------------------------#


fitBeta2 <- function(i,dat){
  dd <- dat[,i]
  #dat1 <- dat2=abs(dd)
  dat1 <- dd[dd>0]
  dat2 <- abs(dd[dd<=0])
  #approximate dist if size is too small
  if(length(dat1)<10)dat1 <- abs(dd)
  if(length(dat2)<10)dat2 <- abs(dd)
  m <- mean(dat1)
  v <- var(dat1)
  shape1_pos <- m * ((m * (1 - m) / v) - 1)
  shape2_pos <- (1 - m) * ((m * (1 - m) / v) - 1)
  m <- mean(dat2)
  v <- var(dat2)
  shape1_neg <- m * ((m * (1 - m) / v) - 1)
  shape2_neg <- (1 - m) * ((m * (1 - m) / v) - 1)
  return(c(shape1_pos,shape2_pos,shape1_neg,shape2_neg))
}

# Handle single-gene sets
handle_single_gene_sets <- function(gSet, gsSize0, es0, rgl, scoreType) {
  idx1 <- which(gsSize0 == 1)
  idx2 <- which(gsSize0 > 1)
  es1 <- es0[idx1]
  es0 <- es0[idx2]
  gsSize0 <- gsSize0[idx2]
  
  null.p <- 1 - (seq_along(rgl) - 1) / (length(rgl) - 1)
  null.n <- -(seq_along(rgl) - 1) / (length(rgl) - 1)
  null.s <- ifelse(null.p > abs(null.n), null.p, null.n)
  null.s.p <- matrix(null.s[null.s > 0], ncol = 1)
  null.s.n <- matrix(null.s[null.s <= 0], ncol = 1)
  
  if (scoreType == "std") {
    p1 <- vapply(es1, function(x) {
      if (x > 0) mean(null.s.p >= x) else mean(null.s.n <= x)
    }, FUN.VALUE = numeric(1))
  } else if (scoreType == "pos") {
    p1 <- vapply(es1, function(x) mean(null.p >= x), FUN.VALUE = numeric(1))
  } else if (scoreType == "neg") {
    p1 <- vapply(es1, function(x) mean(null.n <= x), FUN.VALUE = numeric(1))
  } else {
    stop("scoreType should be std, pos or neg")
  }
  list(p1 = p1, idx1 = idx1, idx2 = idx2, es0 = es0, gsSize0 = gsSize0,
       es00 = c(es1, es0))
}

# Generate null distribution matrix
generate_null_distributions <- function(nSim, gsSize, rgl,ng,ncore,scoreType){
  maxGsSize <- max(gsSize)
  resu <- mclapply(seq_len(nSim), function(k) {
    rgs <- sample.int(ng, maxGsSize)
    vapply(gsSize, function(gsSize_i) {
      rgs.i <- sort(rgs[seq_len(gsSize_i)])
      ngi <- ng - gsSize_i
      x <- rgl[rgs.i]
      sumx <- sum(x)
      y <- rgs.i - seq_len(gsSize_i)
      tmp <- cumsum(x) / sumx - y / ngi
      es.p <- max(tmp)
      es.n <- min(c(0, tmp - x / sum(x)))
      es.s <- ifelse(es.p > abs(es.n), es.p, ifelse(es.p != abs(es.n), es.n, 0))
      switch(scoreType, std = es.s, pos = es.p, neg = es.n)
    }, FUN.VALUE = numeric(1))
  }, mc.cores = ncore)
  do.call(rbind, lapply(resu, unlist))
}

# Calculate permutation p-values
calculate_permutation_pvalues <- function(es0,gsSize0,mat.null,scoreType){
  gsSize0 <- as.character(gsSize0)
  if (scoreType == "pos") {
    vapply(seq_along(es0), function(i) {
      mean(mat.null[, gsSize0[i]] >= es0[i])
    }, FUN.VALUE = numeric(1))
  } else if (scoreType == "neg") {
    vapply(seq_along(es0), function(i) {
      mean(mat.null[, gsSize0[i]] <= es0[i])
    }, FUN.VALUE = numeric(1))
  } else if (scoreType == "std") {
    vapply(seq_along(es0), function(i) {
      null <- mat.null[, gsSize0[i]]
      if (es0[i] > 0) {
        mean(null[null > 0] >= es0[i])
      } else {
        mean(null[null <= 0] <= es0[i])
      }
    }, FUN.VALUE = numeric(1))
  } else {
    stop("scoreType should be pos, neg or std")
  }
}

#Fit beta distribution and adjust p-values
fit_beta_and_adjust_pvalues <- function(es0,mat.null,gsSize0,scoreType,
                                        nSim,permuP){
  if (scoreType == "std") {
    b_para <- do.call(rbind, lapply(seq_len(ncol(mat.null)), fitBeta2,
                                    dat = mat.null))
    rownames(b_para) <- colnames(mat.null)
    b_para <- b_para[as.character(gsSize0), ]
    p <- vapply(seq_along(es0), function(i) {
      if (es0[i] > 0) {
        exp(stats::pbeta(es0[i], b_para[i, 1], b_para[i, 2], log.p = TRUE,
                         lower.tail = FALSE))
      } else {
        exp(stats::pbeta(abs(es0[i]), b_para[i, 3], b_para[i, 4], log.p = TRUE,
                         lower.tail = FALSE))
      }
    }, FUN.VALUE = numeric(1))
    p[p == 0] <- min(p[p != 0])
    p[p == 1] <- max(p[p != 1])
    
    permuP[permuP == 0] <- min(permuP[permuP != 0])
    permuP[permuP == 1] <- max(permuP[permuP != 1])
    
    logP <- -log10(p)
    pfit <- data.frame(p = p, permuP = permuP, es0 = es0)
    pfit <- na.omit(pfit)
    pfit <- pfit[pfit$permuP >= 1 / nSim, ]
    pfit$p <- -log10(pfit$p)
    pfit$permuP <- -log10(pfit$permuP)
    
    model_pos <- lm(permuP ~ p - 1, data = subset(pfit, es0 > 0))
    model_neg <- lm(permuP ~ p - 1, data = subset(pfit, es0 <= 0))
    
    p[es0 > 0] <- predict(model_pos, newdata = data.frame(p = logP[es0 > 0]))
    p[es0 <= 0] <- predict(model_neg, newdata = data.frame(p = logP[es0 <= 0]))
    p <- 1 / 10^p
  } else {
    mat.null <- abs(mat.null)
    b_para <- do.call(rbind, lapply(seq_len(ncol(mat.null)), fitBeta,
                                    dat = mat.null))
    rownames(b_para) <- colnames(mat.null)
    b_para <- b_para[as.character(gsSize0), ]
    es0 <- abs(es0)
    p <- vapply(seq_along(es0), function(i) {
      exp(stats::pbeta(es0[i], b_para[i, 1], b_para[i, 2], log.p = TRUE, 
                       lower.tail = FALSE))
    }, FUN.VALUE = numeric(1))
    p[p == 0] <- min(p[p != 0])
  }
  p
}


# betaFitPermu function
betaFitPermu <- function(gSet, rgl, nSim = 10000, es0, ncore, scoreType) {
  gsSize0 <- vapply(gSet, length, FUN.VALUE = numeric(1))
  ng <- length(rgl)
  
  if (1 %in% gsSize0) {
    result <- handle_single_gene_sets(gSet, gsSize0, es0, rgl, scoreType)
    p1 <- result$p1
    idx1 <- result$idx1
    idx2 <- result$idx2
    es0 <- result$es0
    gsSize0 <- result$gsSize0
    es00 <- result$es00
  }
  
  gsSize <- unique(gsSize0)
  mat.null <- generate_null_distributions(nSim,gsSize,rgl,ng,ncore,scoreType)
  colnames(mat.null) <- gsSize
  
  permuP <- calculate_permutation_pvalues(es0, gsSize0, mat.null, scoreType)
  p <- fit_beta_and_adjust_pvalues(es0, mat.null, gsSize0, scoreType,
                                   nSim,permuP)
  
  if (exists("es00")) {
    p2 <- p; p <- rep(NA, length(es00)); p[idx1] <- p1; p[idx2] <- p2
    p2 <- permuP; permuP <- rep(NA, length(es00)); permuP[idx1] <- p1; 
    permuP[idx2] <- p2
  }
  
  data.frame(p = p, permuP = permuP)
}

#main function
gsRank <- function(stats,outfile="gsRank",scoreType="std",gSetName=NULL,
                   geneSet=NULL,gseaParam=1,species="Human",nperm=1e4,ncore=1){
  
  genep <- stats
  genep <- genep[!is.na(genep)]
  genep <- sort(genep,decreasing=TRUE)
  genep <- abs(genep)
  genep <- genep^gseaParam
  
  gSetName12 <- ParseGsetName(geneSet,gSetName)
  gSetName1 <- gSetName12$gSetName1
  gSetName2 <- gSetName12$gSetName2
  
  for(gset1 in gSetName1){
    for(gset2 in gSetName2){
      
      if(is.null(geneSet)){
        if(gset1=="KEGG"){geneSet <- getKEGG(species)
        }else if(gset1=="GO"){geneSet <- getGO(gset2,species)
        }else if(gset1=="MSigDB"){geneSet <- getMSigDB(gset2,species)
        }else if(gset1=="Reactome"){geneSet <- getReactome(species)
        }else{stop("Please specify gSetName as KEGG, GO, MSigDB, 
            or Reactome, or provide a geneSet list")}
      }
      
      # Check geneSet, remove genes that are not in test universe
      universe <- names(genep)
      if(!is.list(geneSet))stop("geneSet should be a list object")
      geneSet <- lapply(geneSet,function(x){x <- as.character(as.vector(x)) 
      x <- x[!is.na(x)];x <- x[x %in% universe];unique(x)})
      inUniv <- vapply(geneSet, function(x) length(x) > 0,FUN.VALUE=logical(1))
      geneSet <- geneSet[inUniv]
      
      es0 <- vapply(geneSet,function(gsi){
        es1 <- es(gsi,genep)
        switch(scoreType,"std"={
          m1 <- max(es1); m2 <- min(es1)
          ifelse(m1 > abs(m2), m1, ifelse(m1 == abs(m2), 0, m2))
        },"pos"=max(es1),"neg"=min(es1))},FUN.VALUE=numeric(1))
      p <- betaFitPermu(gSet=geneSet,rgl=genep,nSim=nperm,es0=es0,ncore=ncore,
                        scoreType=scoreType)
      nGene <- vapply(geneSet,length,FUN.VALUE=numeric(1))
      resu <- data.frame(ID=names(geneSet),nGene=nGene,EnrichScore=es0,p=p$p,
                         permuP=p$permuP)
      resu$FDR <- p.adjust(resu$p,method="fdr")
      
      outfn <- paste0(outfile,"_",gset1,"_",gset2,".csv")
      write.csv(resu,outfn,row.names=FALSE)
    }}}

