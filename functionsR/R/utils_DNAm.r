#' Extract scan date from idat file
#' 
#'
#' @param idat file path to a single idat file from either colour channel for a sample e.g. 205096370026_R02C01_Red.idat
#' @return date idat file was scanned by 
#'
#' @importFrom illuminaio readIDAT
#' @author E Hannon
#' @export
getScanDate<-function(idat){
	tmp<-readIDAT(idat)$RunInfo
	return(unique(unlist(lapply(strsplit(tmp[which(tmp[,"BlockType"] == "Scan"),1], " "), head, n = 1)))[1])
}

#' Reformat beta matrix so rownames ordered to match a specificed list
#' 
#'
#' @param bMat matrix of beta values
#' @param pList a vector of rownames in the required order
#' @return a reordered matrix of beta values
#'
#' @author E Hannon
#' @export
reformatBetas <- function(bMat, pList){
	index<-match(pList, rownames(bMat))
	return(bMat[index,])
}


#' Calculate lambda from genome-wide association results
#' 
#'
#' @param pvals a vector of p-values to estimate lambda from
#' @return the lambda value
#'
#' @author E Hannon
#' @export
estlambda<-function(pvals){
	z = qnorm(pvals / 2)
	## calculates lambda
	lambda = round(median(z^2) / 0.454, 3)
	return(lambda)
}