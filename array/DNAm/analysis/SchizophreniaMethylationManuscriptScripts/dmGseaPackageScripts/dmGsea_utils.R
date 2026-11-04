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



#----------------------------------------------------------------------#
# C+P Functions
#----------------------------------------------------------------------#


#es: Enrichment Score;gs: a gene set; rgl:ranked gene test result
es <- function(gsi,rgl){
  ng <- length(rgl)
  idx <- names(rgl) %in% gsi
  x <- rep(0,ng);x[idx] <- 1;y <- 1-x
  cumsum(x*rgl/sum(rgl[idx])-y/(ng-length(gsi)))
}

fitBeta <- function(i,dat){
  m <- mean(dat[,i])
  v <- var(dat[,i])
  shape1 <- m * ((m * (1 - m) / v) - 1)
  shape2 <- (1 - m) * ((m * (1 - m) / v) - 1)
  return(c(shape1,shape2))
}

ParseGsetName <- function(geneSet=NULL,gSetName=NULL){
  if(is.null(geneSet)){
    if(is.null(gSetName)){stop("gSetName need to be specified")}else{
      tmp <- unlist(strsplit(gSetName,"\\|"))
      if(length(tmp)==1){
        if(tmp[1] %in% c("KEGG","Reactome")){gSetName1 <- tmp[1];gSetName2<- tmp[1]
        }else if(tmp[1]=="GO"){gSetName1 <- tmp[1];gSetName2 <- c("BP","MF","CC")
        }else if(tmp[1]=="MSigDB"){gSetName1 <- tmp[1];gSetName2 <- c("H","C1",
                                                                      "C2","C3","C4","C5","C6","C7","C8")
        }else{stop("We can not find ",gSetName," currently, please specify the
        pathway list using geneSet" )
        }
        
      }else{
        gSetName1 <- tmp[1];gSetName2 <- tmp[2]
      }
    }
  }else{gSetName1 <- "userSet";gSetName2 <- "userSet"}
  return(list(gSetName1=gSetName1,gSetName2=gSetName2))
}



getGO <- function(subset="BP",species="Human"){
  if(species=="Human"){
    if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("org.Hs.eg.db package required but not installed.")
    if(!requireNamespace("GO.db", quietly = TRUE))
      stop("GO.db package required but not installed.")
    egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", 
                                           "org.Hs.eg.db")
  }else if(species=="Mouse"){
    if(!requireNamespace("org.Mm.eg.db", quietly = TRUE))
      stop("org.Mm.eg.db package required but not installed.")
    egGO2ALLEGS <- utils::getFromNamespace("org.Mm.egGO2ALLEGS", 
                                           "org.Mm.eg.db")
  }else{stop("We only have Human and Mouse data currently, please
            provide geneSet list")}
  GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", 
                                                          "go_id", "Ontology")]
  GeneID.PathID <- GeneID.PathID[GeneID.PathID$Ontology==subset,]
  GeneID.PathID <- unique(GeneID.PathID)
  go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
  GOID.TERM <- AnnotationDbi::select(GO.db::GO.db, 
                                     keys=unique(GeneID.PathID$go_id), 
                                     columns=c("GOID","ONTOLOGY","TERM"), 
                                     keytype="GOID")
  GOID.TERM <- GOID.TERM[match(names(go),GOID.TERM[,1]),]
  names(go) <- paste(GOID.TERM$GOID,GOID.TERM$TERM,sep="|")
  return(go)
}

getKEGG <- function(species="Human"){
  if(!requireNamespace("KEGGREST",quietly=TRUE))stop("KEGGREST package
                                required")
  if(species=="Human"){species <- "hsa"
  }else if(species=="Mouse"){species <- "mmu"
  }else{stop("We only have Human and Mouse data currently, please provide
            geneSet list")}
  
  links    <- KEGGREST::keggLink("pathway", species)
  path <- data.frame(geneid=names(links),pathwayid=links)
  path$geneid <- sub(".*:","",path$geneid)
  path$pathwayid <- sub(".*:","",path$pathwayid)
  kegg <- tapply(path$geneid, path$pathwayid,list)
  name <- KEGGREST::keggList("pathway", species)
  names(kegg) <- paste0(names(kegg),"|",name[names(kegg)])
  return(kegg)
}

getMSigDB <- function(subset="C2",species="Human"){
  if(species=="Human"){species <- "Homo sapiens"
  }else if(species=="Mouse"){species <- "Mus musculus"
  }else{"We only have Human and Mouse data currently, please provide 
    geneSet list"}
  if(!requireNamespace("msigdbr",quietly=TRUE))stop("msigdbr package 
                                required but not installed")
  gsets <- msigdbr::msigdbr(species = species, category = subset)
  gsets$id <- paste(gsets$gs_id,gsets$gs_subcat,gsets$gs_name,sep="|")
  gsets <- tapply(gsets$entrez_gene,gsets$id,list)
  return(gsets)
}

getReactome <- function(species="Human"){
  if(species=="Human"){species <- "Homo sapiens"
  }else if(species=="Mouse"){species <- "Mus musculus"
  }else{"We only have Human and Mouse data currently, please provide
    geneSet list"}
  
  #    URL <- "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt"
  #    n2r <- read.table(URL, sep = "\t", quote = "\"", fill = TRUE,
  #                comment.char = "", stringsAsFactors = FALSE)
  n2r <- readRDS(system.file("ReactomeDB.rds",package="dmGsea"))
  n2r <- n2r[n2r$V6 == species,]
  n2r$id <- paste(n2r$V2,n2r$V4,sep="|")
  gsets <- tapply(n2r$V1,n2r$id,list)
  return(gsets)
}


getIlluminaAnnotation <- function(arrayType=c("450K","EPIC"))
{
  if(arrayType=="450K"){
    annopkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  }else if(arrayType=="EPIC") {
    annopkg <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  }else {stop("arrayType should be 450K or EPIC ")}
  
  if(!requireNamespace(eval(annopkg),quietly=TRUE))stop(annopkg,
                                                        " required but not installed")
  anno <- minfi::getAnnotation(eval(annopkg))
  
  # get rid of the non-CpG sites and CpGs that are not annotated
  ann.keep<- subset(anno,grepl("^cg",anno$Name) & anno$UCSC_RefGene_Name!="")
  # get individual isoforms for each CpG
  isolist <- strsplit(ann.keep$UCSC_RefGene_Accession,split=";")
  grouplist <- strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  index <- rep(seq_len(nrow(ann.keep)),times=vapply(isolist,length,
                                                    FUN.VALUE=numeric(1)))
  GeneProbeTable <- data.frame(iso=unlist(isolist),group=unlist(grouplist),
                               ann.keep[index,c("Name","chr","pos")],stringsAsFactors=FALSE)
  if(!requireNamespace("org.Hs.eg.db",quietly=TRUE))stop(
    "org.Hs.eg.db package required")
  eg <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
                              keys=GeneProbeTable$iso, 
                              columns="ENTREZID", 
                              keytype="ACCNUM")
  
  GeneProbeTable$entrezid <- eg$ENTREZID
  GeneProbeTable <- GeneProbeTable[!is.na(GeneProbeTable$entrezid),]
  # keep unique cpg by gene id
  id <- paste(GeneProbeTable$Name,GeneProbeTable$entrezid,sep=".")
  d <- duplicated(id)
  GeneProbeTable <- GeneProbeTable[!d,]
  return(GeneProbeTable)
}


# gsea test 
p.gs <- function(gs,gene.info){
  gs.ng <- length(gs)
  gs.score <- sum(gene.info$score[gene.info$entrezid %in% gs])
  #    sig.gene <- paste(gene.info$entrezid[gene.info$entrezid %in% gs &
  # gene.info$score==1],collapse=";")
  n.sg <- sum(gene.info$score)
  if(gs.score==0){p <- 1}else{    
    pw.red <- with(gene.info,mean(pwf[entrezid %in% gs]))
    pw.white <- with(gene.info,mean(pwf[!(entrezid %in% gs)]))
    #    odds <- pw.red*(1-pw.white)/(pw.white*(1-pw.red))
    #    odds <- pw.red/pw.white
    odds <- min((1-pw.red)*pw.white/((1-pw.white)*pw.red),1e+100)
    
    ng2 <- nrow(gene.info)-gs.ng
    p <- BiasedUrn::pFNCHypergeo(gs.score,gs.ng,ng2,n.sg,odds,lower.tail=FALSE)
    +BiasedUrn::dFNCHypergeo(gs.score,gs.ng,ng2,n.sg,odds)
  }
  c(gs.ng,gs.score,p)
}

# #combine P values for each gene or gene group
# combinep <- function(gid,Data4Cor,GeneProbeTablep,combpMethod,combpAdjust){
#   #    if(!requireNamespace("poolr",quietly=TRUE))stop("poolr package required")
#   subset <- GeneProbeTablep[GeneProbeTablep$entrezid==gid,]
#   allp <- subset$p
#   if(length(allp)>=2){
#     if(is.null(Data4Cor)){
#       eval(parse(text=paste0("poolr::",combpMethod,"(p=allp,adjust='none')")))$p
#     }else{
#       dat <- t(Data4Cor[as.vector(subset$Name),])
#       Rmatrix <- cor(dat,use="pairwise.complete.obs")
#       diag(Rmatrix) <- rep(1, nrow(Rmatrix))
#       eval(parse(text=paste0("poolr::",combpMethod,"(p=allp,adjust=combpAdjust,
#             R=Rmatrix)")))$p
#     }}else{allp}
# }


geneID2geneName    <- function(gid=NULL,species="Human"){
  if(species=="Human"){
    gene_names <- mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gid,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  }else if(species=="Mouse"){
    gene_names <- mapIds(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = gid,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  }else{
    gene_names <- gid
    names(gene_names) <- gid
  }
  gene_names[is.na(gene_names)] <-
    paste0("geneid:",names(gene_names)[is.na(gene_names)])
  return(gene_names)
}




#----------------------------------------------------------------------#