
## load config variables
if (is.na(params$project)){
  project<-"rizzardi"
  source("ATACSeq/config/config.r")

} else {
    project<-params$project
  source(params$configFile)
}
