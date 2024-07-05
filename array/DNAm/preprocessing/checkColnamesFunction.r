

# Checks sample sheet columns are formatted correctly
checkColnames <- function(sampleSheet, exp_cols, type, verbose=T){
	metrics <- list()
	
	all_pres <- all(exp_cols %in% colnames(sampleSheet))
	metrics[['allPresent']] <- all_pres
	
	if(all_pres){
		if(verbose){print(paste(type, "column names: All present"))}
	}else{
	
		if(verbose){print(paste(type, "column names: Not all present"))}
		
		#a) find correct cols
		pres <- exp_cols[exp_cols %in% colnames(sampleSheet)]
		if(length(pres)>0){
			if(verbose){print(c("Identical colnames:",pres))}
		}
		metrics[['presentIDs']] <- pres
		
		#b) find absent cols
		ab <- exp_cols[exp_cols %ni% colnames(sampleSheet)]
		if(verbose){print(c("Absent colnames:",ab))}
		metrics[['absentIDs']] <- ab
		
		#c) find close matches
		near <- exp_cols[amatch(colnames(sampleSheet), exp_cols, maxDist = 2)] # maxDist = number of characters difference

		if(length(near[-which(is.na(near))])>0){
			nrmatch <- data.frame(Expected.colname=near, Received.colname=colnames(sampleSheet))
			nrmatch <- nrmatch[which(!is.na(nrmatch$Expected)),]
			# when near==expected, remove this match and any close matches also to this colname
			equal <- which(as.character(nrmatch$Expected) == as.character(nrmatch$Received.colname))
			if(length(equal)>0){
				other <- which(as.character(nrmatch$Expected) %in% as.character(nrmatch$Expected)[equal]) # near matches to a colname that is present, also captures equal
				nrmatch <- nrmatch[-other,]
			}
			if(nrow(nrmatch)>0){
				if(verbose){print(paste("Identified", nrow(nrmatch), "near match(es)..."))}
				if(verbose){print(nrmatch)}
				metrics[['nearMatchIDs']] <- nrmatch
			}
		}
		
	}
	if(!verbose){
		return(metrics)
	}
}