# editted version of enrichPeakOverlap from CHIPSeeker package to 
# 1) report FC statistic
# 2) test both over and under enrichment.
# calculate fold change statistic as ratio of obeserved overlap to mean overlap.

#' calculate overlap significant of ChIP experiments based on the genome coordinations
#'
#'
#' @title enrichPeakOverlap
#' @param queryPeak query bed file or GRanges object
#' @param targetPeak target bed file(s) or folder that containing bed files or a list of GRanges objects
#' @param TxDb TxDb
#' @param pAdjustMethod pvalue adjustment method
#' @param nShuffle shuffle numbers
#' @param chainFile chain file for liftOver
#' @param pool logical, whether pool target peaks
#' @param mc.cores number of cores, see \link[parallel]{mclapply}
#' @param verbose logical
#' @return data.frame
#' @export
#' @importFrom rtracklayer import.chain
#' @importFrom rtracklayer liftOver
#' @author G Yu
enrichPeakOverlap <- function(queryPeak, targetPeak, TxDb=NULL, pAdjustMethod="BH", nShuffle=1000,
                              chainFile=NULL, pool=TRUE, mc.cores=detectCores()-1, verbose=TRUE) {
    TxDb <- loadTxDb(TxDb)
    query.gr <- loadPeak(queryPeak)
    if (is(targetPeak[1], "GRanges") || is(targetPeak[[1]], "GRanges")) {
        target.gr <- targetPeak
        targetFiles <- NULL
    } else {
        targetFiles <- parse_targetPeak_Param(targetPeak)
        target.gr <- lapply(targetFiles, loadPeak)
    }

    if (!is.null(chainFile)) {
        chain <- import.chain(chainFile)
        target.gr <- lapply(target.gr, liftOver, chain=chain)
    }

    if (pool) {
        p.ol <- enrichOverlap.peak.internal(query.gr, target.gr, TxDb, nShuffle,
                                            mc.cores=mc.cores,verbose=verbose)
    } else {
        res_list <- lapply(1:length(target.gr), function(i) {
            enrichPeakOverlap(queryPeak = queryPeak,
                              targetPeak = target.gr[i],
                              TxDb = TxDb,
                              pAdjustMethod = pAdjustMethod,
                              nShuffle = nShuffle,
                              chainFile = chainFile,
                              mc.cores = mc.cores,
                              verbose = verbose)
        })
        res <- do.call("rbind", res_list)
        return(res)
    }

    if (is.null(p.ol$pvalueOver)) {
        p2 <- p1 <- padj <- NA
    } else {
        p1 <- p.ol$pvalueOver
		p2 <- p.ol$pvalueUnder
        padj <- p.adjust(p1, method=pAdjustMethod)
    }

    ol <- p.ol$overlap
	fc <- p.ol$foldchange
	fc.l95 <- p.ol$foldchangeL95
	fc.u95 <- p.ol$foldchangeU95

    if (is(queryPeak, "GRanges")) {
        qSample <- "queryPeak"
    } else {
        # remove path, only keep file name
        qSample <- basename(queryPeak)
    }

    if (is.null(targetFiles)) {
        tSample <- names(target.gr)
        if(is.null(tSample)) {
            tSample <- paste0("targetPeak", seq_along(target.gr))
        }
    } else {
        tSample <- basename(targetFiles)
    }

    res <- data.frame(qSample=qSample,
                      tSample=tSample,
                      qLen=length(query.gr),
                      tLen=unlist(lapply(target.gr, length)),
                      N_OL=ol,
					  foldchange=fc,
					  foldchangeL95=fc.l95,
					  foldchangeU95=fc.u95,
                      pvalueOver=p1,
					  pvalueUnder=p2,
                      p.adjust=padj)

    return(res)
}

#' @import GenomeInfoDb
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
enrichOverlap.peak.internal <- function(query.gr, target.gr, TxDb, nShuffle=1000, mc.cores=detectCores()-1, verbose=TRUE) {
    if (verbose) {
        cat(">> permutation test of peak overlap...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    }

    idx <- sample(1:length(target.gr), nShuffle, replace=TRUE)
    len <- unlist(lapply(target.gr, length))

    if(Sys.info()[1] == "Windows") {
        qLen <- lapply(target.gr, function(tt) {
            length(intersect(query.gr, tt))
        })
    } else {
        qLen <- mclapply(target.gr, function(tt) {
            length(intersect(query.gr, tt))
        }, mc.cores=mc.cores)
    }
    qLen <- unlist(qLen)
    # query ratio
    qr <- qLen/len
	

    if (nShuffle < 1) {
        res <- list(pvalueOver=NULL, pvalueUnder=NULL, overlap=qLen, foldchange=NULL, foldchangeL95 = NULL, foldchangeU95 = NULL)
        return(res)
    }

    if (verbose) {
        pb <- txtProgressBar(min=0, max=nShuffle, style=3)
    }
    if(Sys.info()[1] == "Windows") {
        rr <- lapply(seq_along(idx), function(j) {
            if (verbose) {
                setTxtProgressBar(pb, j)
            }
            i <- idx[j]
            tarShuffle <- shuffle(target.gr[[i]], TxDb)
            length(intersect(query.gr, tarShuffle))/len[i]
        })
    } else {
        rr <- mclapply(seq_along(idx), function(j) {
            if (verbose) {
                setTxtProgressBar(pb, j)
            }
            i <- idx[j]
            tarShuffle <- shuffle(target.gr[[i]], TxDb)
            length(intersect(query.gr, tarShuffle))/len[i]
        }, mc.cores=mc.cores,mc.preschedule = FALSE)
    }

    if (verbose) {
        close(pb)
    }

    rr <- unlist(rr) # random ratio
	
	fc<-unlist(lapply(qr, function(q) q/mean(rr)))
	fc.l95<-unlist(lapply(qr, function(q) q/quantile(rr, 0.975)))
	fc.u95<-unlist(lapply(qr, function(q) q/quantile(rr, 0.025)))
	
    # p <- lapply(qr, function(q) mean(rr>q))
	# this is calculation of one-side p for over enrichment
    p1 <- unlist(lapply(qr, function(q) (sum(rr>q)+1)/(length(rr)+1)))
	# add in calculation of one-sided p for under enrichment
	p2 <- unlist(lapply(qr, function(q) (sum(rr<q)+1)/(length(rr)+1)))
	
    res <- list(pvalueOver=p1, pvalueUnder=p2, overlap=qLen, foldchange=fc, foldchangeL95 = fc.l95, foldchangeU95 = fc.u95)
    return(res)
}


# additional utility functions that need to be loaded but not editted

loadPeak <- function(peak, verbose=FALSE) {
    if (is(peak, "GRanges")) {
        peak.gr <- peak
    } else if (file.exists(peak)) {
        if (verbose)
            cat(">> loading peak file...\t\t\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        peak.gr <- readPeakFile(peak, as="GRanges")
    } else {
        stop("peak should be GRanges object or a peak file...")
    }
    return(peak.gr)
}

#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
loadTxDb <- function(TxDb) {
    if ( is.null(TxDb) ) {
        warning(">> TxDb is not specified, use 'TxDb.Hsapiens.UCSC.hg19.knownGene' by default...")
        TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    return(TxDb)
}

# function to count min overlap
countOverlapMin<-function(query, test, min = 0.5){
	# if provided as proportion then looking for overlap proportional to query region
	if(min < 1){
		whichOverlap<-findOverlaps(test,query)
		testLengths<-width(test[queryHits(whichOverlap)])
		thresLength<-min*testLengths
		
	} else {
		thresLength = min
	}
	return(sum(width(intersect(test, query)) > thresLength))
}