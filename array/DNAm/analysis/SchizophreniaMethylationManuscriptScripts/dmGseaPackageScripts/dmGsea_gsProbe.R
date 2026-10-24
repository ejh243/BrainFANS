# Get GSEA CpG data
get.gsea.cpg <- function(cpg.p,cpg.gene,nTopProbe=NULL,FDRthre=0.05){
  all.cpg <- cpg.p$Name; n.cpg <- nrow(cpg.p)
  
  #significant CpGs
  if(!is.null(nTopProbe)){sig.cpg <- cpg.p$Name[order(cpg.p$p,
                                                      decreasing = FALSE)][seq_len(nTopProbe)]}else{
                                                        sig.cpg <- with(cpg.p,Name[p.adjust(p,method="fdr")<FDRthre])}
  
  n.sig.cpg <- length(sig.cpg)
  all.gene <- unique(cpg.gene$entrezid)
  sig.gene <- with(cpg.gene,unique(entrezid[Name %in% sig.cpg]))
  
  if(length(sig.gene)<10){stop("Number of significant gene is less than 10,
                please adjust threshold or use ranking method in gsGene()")}
  # pwf from cpg
  gene.freq <- c(table(cpg.gene$entrezid)[all.gene])
  gene.prob <- dhyper(0,gene.freq,n.cpg-gene.freq,n.sig.cpg)
  gene.prob[gene.prob==0] <- min(gene.prob[!gene.prob==0])
  
  data.frame(entrezid=all.gene,pwf=unname(gene.prob),score=ifelse(
    all.gene %in% sig.gene,1,0))
}

# Validate and prepare input data
prep.gsProbe.input <- function(probe.p, GeneProbeTable, arrayType, sigProbe,
                               allProbe, nTopProbe) {
  if (is.null(GeneProbeTable)) {
    if (!(arrayType %in% c("450K", "EPIC")))
      stop("Specify arrayType as 450K or EPIC, or provide GeneProbeTable")
    GeneProbeTable <- getIlluminaAnnotation(arrayType = arrayType)
  }
  if (!is.null(sigProbe)) {
    if (is.null(allProbe)) allProbe <- unique(GeneProbeTable$Name)
    probe.p <- data.frame(Name = allProbe, p = 0.999)
    probe.p$p[probe.p$Name %in% sigProbe] <- 1e-10
    nTopProbe <- sum(probe.p$p < 0.05)
  }
  names(probe.p)[which(names(probe.p) == "probe")] <- "Name"
  probe.p$Name <- as.character(as.vector(probe.p$Name))
  GeneProbeTable$p <- probe.p$p[match(GeneProbeTable$Name, probe.p$Name)]
  GeneProbeTable <- GeneProbeTable[complete.cases(GeneProbeTable), ]
  GeneProbeTable <- GeneProbeTable[order(GeneProbeTable$entrezid,
                                         GeneProbeTable$p), ]
  GeneProbeTable$Name <- as.character(as.vector(GeneProbeTable$Name))
  GeneProbeTable$entrezid <- as.character(as.vector(GeneProbeTable$entrezid))
  GeneProbeTable$p <- as.numeric(as.vector(GeneProbeTable$p))
  list(probe.p = probe.p, GeneProbeTable=GeneProbeTable,nTopProbe=nTopProbe)
}

# Run GSEA and write results
run.gsea.probe <- function(probe.p, GeneProbeTable, nTopProbe, FDRthre, geneSet,
                           gSetName12, ncore, species, outfile) {
  cpg.p <- GeneProbeTable[!duplicated(GeneProbeTable$Name), ]
  cpg.gene <- GeneProbeTable[, c("Name", "entrezid")]
  gsea.cpg <- get.gsea.cpg(cpg.p, cpg.gene, nTopProbe, FDRthre = FDRthre)
  for (gset1 in gSetName12$gSetName1) {
    for (gset2 in gSetName12$gSetName2) {
      if (is.null(geneSet)) {
        if (gset1 == "KEGG") geneSet <- getKEGG(species)
        else if (gset1 == "GO") geneSet <- getGO(gset2, species)
        else if (gset1 == "MSigDB") geneSet <- getMSigDB(gset2, species)
        else if (gset1 == "Reactome") geneSet <- getReactome(species)
        else stop("Specify gSetName as KEGG, GO, MSigDB, or Reactome")
      }
      universe <- unique(GeneProbeTable$entrezid)
      if (!is.list(geneSet)) stop("geneSet should be a list object")
      geneSet <- lapply(geneSet, function(x) {
        x <- as.character(as.vector(x))
        x <- x[!is.na(x)]
        x <- x[x %in% universe]
        unique(x)
      })
      inUniv <- vapply(geneSet, function(x) length(x) > 0,
                       FUN.VALUE = logical(1))
      geneSet <- geneSet[inUniv]
      resu1 <- mclapply(geneSet, function(x) p.gs(x, gsea.cpg),
                        mc.cores = ncore)
      resu <- do.call(rbind, lapply(resu1, unlist))
      resu <- data.frame(ID = names(geneSet), resu)
      names(resu) <- c("ID", "nGene", "nSigGene", "p")
      sig.gene <- unique(gsea.cpg$entrezid[gsea.cpg$score == 1])
      sgName <- geneID2geneName(sig.gene, species)
      resu$Sig.Gene <- vapply(geneSet, function(x)
        paste(sgName[sig.gene[sig.gene %in% x]], collapse = ";"),
        FUN.VALUE = character(1))
      resu$FDR <- p.adjust(resu$p, method = "fdr")
      outfn <- paste0(outfile, "_", gset1, "_", gset2, ".csv")
      write.csv(resu, outfn, row.names = FALSE)
    }
  }
}

# Main function
gsProbe <- function(probe.p, FDRthre = 0.05, nTopProbe = NULL, sigProbe = NULL,
                    allProbe = NULL, GeneProbeTable = NULL, arrayType = NULL,
                    gSetName = NULL, geneSet = NULL, species = "Human",
                    outfile = "gsProbe", ncore=1) {
  input <- prep.gsProbe.input(probe.p, GeneProbeTable, arrayType, sigProbe,
                              allProbe, nTopProbe)
  probe.p <- input$probe.p
  GeneProbeTable <- input$GeneProbeTable
  nTopProbe <- input$nTopProbe
  gSetName12 <- ParseGsetName(geneSet, gSetName)
  run.gsea.probe(probe.p, GeneProbeTable, nTopProbe, FDRthre, geneSet,
                 gSetName12, ncore, species, outfile)
}