## Written by Eilis
## Design Methylation array plate layout profiling
## 1 case 1 control per chip that were facs sorted on the same run/day. 
## randomise the order of sample type

setwd("")
sampleSheet<-read.table("MRC2Chip loc_.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
chip.positions<-paste("R0", 1:8, "C01", sep = "")

## start with dates where 2 IDs were processed
outOrder<-NULL
dates<-sample(names(which(table(sampleSheet$Sort.Date) > 6)))
for(each in dates){
  IDs<-unique(sampleSheet$Individual[which(sampleSheet$Sort.Date == each)])
  # mix up  cases and controls
  IDs<-sample(IDs)
  for(entry in IDs){
    selected<-which(sampleSheet$Individual == entry)
    ## if sample missing replace with NA for now
    while(length(selected) < 4){
        selected<-c(selected, NA)
    }
    selected<-sample(selected)
    outOrder<-c(outOrder, selected)
  }
}

## add on samples run on their own, as number of cases and controls equal we can run one of each per chip

dates<-sample(names(which(table(sampleSheet$Sort.Date) < 5)))
datesByPheno<-unique(sampleSheet[which(sampleSheet$Sort.Date %in% dates), c("Sort.Date", "Phenotype")])
dates.con<-sample(datesByPheno$Sort.Date[which(datesByPheno$Phenotype == "Control")])
dates.case<-sample(datesByPheno$Sort.Date[which(datesByPheno$Phenotype == "Schizophrenia")])

## alternate order of cases and controls
caseConOrder<-c(1,2)

for(i in 1:length(dates.con)){
  caseConSelect<-sample(caseConOrder,1)
  if(caseConSelect == 1){
    IDs<-unique(sampleSheet$Individual[which(sampleSheet$Sort.Date == dates.con[i])])
    selected<-which(sampleSheet$Individual == IDs[1])
    ## if sample missing replace with NA for now
    while(length(selected) < 4){
       selected<-c(selected, NA)
    }
    selected<-sample(selected)
    outOrder<-c(outOrder, selected)
    }

  IDs<-unique(sampleSheet$Individual[which(sampleSheet$Sort.Date == dates.case[i])])
  selected<-which(sampleSheet$Individual == IDs[1])
  ## if sample missing replace with NA for now
  while(length(selected) < 4){
     selected<-c(selected, NA)
  }
  selected<-sample(selected)
  outOrder<-c(outOrder, selected)
  
  if(caseConSelect == 2){
    IDs<-unique(sampleSheet$Individual[which(sampleSheet$Sort.Date == dates.con[i])])
    
    selected<-which(sampleSheet$Individual == IDs[1])
    ## if sample missing replace with NA for now
    while(length(selected) < 4){
      selected<-c(selected, NA)
    }
    selected<-sample(selected)
    outOrder<-c(outOrder, selected)
  }
}

## check all samples included
if(length(which(!1:nrow(sampleSheet) %in% outOrder)) > 0){
  print("ERROR: Some samples not included")
}

sampleSheet<-sampleSheet[outOrder,]
sampleSheet$Chip.Location<-chip.positions

write.table(sampleSheet, "SCZDNAmethFACsSortedPlateDesign14052019.txt", row.names = FALSE)

