
## load config variables
if (is.na(params$project)){
  project<-"rizzardi"
  source("ATACSeq/config/config.r")

} else {
    project<-params$project
  source(params$configFile)
}


## create colourblind friendly palette
colorBlindGrey8   <- c("#009E73", "#CC79A7", "#D55E00", "#999999", 
                       "#F0E442", "#0072B2",  "#E69F00", "#56B4E9")

