clump<-function(dat, window = 5000, thres1 = 1e-7, thres2 = 5e-5, listProbes = TRUE){
#' takes all significant association less than thres1 in order of significance as an index site and identifies all other associations within specified window less than thres2. If a less significant probe is clumped with a more significant probe it is excluded for consideration as an index site. 
#'
#' @param dat A matrix with the columns ids, pval, chr and pos,
#' @param window Distance in bp of region around site to look for other sinificant associations.
#' @param thres1 P value threshold to select significant sites to clump.
#' @param thres2 P value threshold to select additional sites to clump with site under consideration
#' @return a ma 

	tmp<-dat[order(dat$pval),]
	out<-NULL
	while(tmp[1,2] < thres1){
		chrom<-tmp[1,3]
		start<-tmp[1,4]-window
		end<-tmp[1,4]+window
		index<-which(tmp$chr == chrom & tmp$pos <= end & tmp$pos >= start & tmp$pos != tmp[1,4])	## remove self from consideration

		n.sig<-length(which(tmp[index,2] < thres2))
		if(n.sig > 0){
			probes<-tmp[index,1][which(tmp[index,2] < thres2)]
			if(n.sig > 1){
				probes<-paste(probes, collapse = ";")
				}
		} else {
			probes<-NA
		}
		out<-rbind(out, c(unlist(tmp[1,]), length(index), n.sig, probes))
		## remove these probes from furture consideration
		tmp<-tmp[-c(1,index),]
	}
	colnames(out)<-c("Index", "Pvalue", "Chr", "Position", "nSitesinWindow", "nSigSites", "SigDNAmSites")
	return(out)
	
}
