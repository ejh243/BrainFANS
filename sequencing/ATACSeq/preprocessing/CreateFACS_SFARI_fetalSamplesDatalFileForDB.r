#Emma Walker
#19/01/2022
#E.M.Walker@exeter.ac.uk


# set up
setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/EmmaW/ATAC_misc/")

# read in lab sheet from Jon Davies
samps <- read.csv("FACSATAC_SampleSheet_10518_10587_EW.csv", stringsAsFactors = F)

# remove adult samples
samps <- samps[-c(which(samps$Adult.Fetal == "Adult")),]

#cols that should contain identical info within an individual
cols_same <- c("Individual.ID", "Adult.Fetal", "Age", "Sex", "Institute", "Brain.Region", "Ethnicity")


#check if any samples have different info for phenotype data that should be the same (e.g. age)

IDs <- unique(samps$Individual.ID)

for(i in 1:length(IDs)){
  x <- samps[which(samps$Individual.ID == IDs[i]),]
  for(j in cols_same){
    if(length(unique(x[,j])) > 1){
      print(paste0("FAIL: ", j, " sample info is not the same for sample ", x$Individual.ID[1]))
    }
  }}

# test code above works by simulating an 'x' object with multiple pheno data in Adult.Fetal col
#i = 1
#x <- samps[which(samps$Individual.ID == IDs[i]),]
#x$Adult.Fetal[3] <- "MISS"
#for(j in cols_same){
  #if(length(unique(x[,j])) > 1){
    #print(paste0("FAIL: ", j, " sample info is not the same for sample ", x$Individual.ID[1]))
  #}
#}

# 3 samples fail the check: 
#[1] "FAIL: Brain.Region sample info is not the same for sample 11834"
#[1] "FAIL: Brain.Region sample info is not the same for sample 13781"
#[1] "FAIL: Ethnicity sample info is not the same for sample 13360"

# manually add in ethnicity for 13360_SATB2neg as 'SAS' to match 13360_SATB2pos in original csv from JonD
# add note in note col for differint brain regions in 11834 and 13781

# note Alice has an alternative ID for NP12/078 (alt_ID R12/078) was also used in the July 2021 cohort of fetal FACS DNAm


# fill in column info for samples
dat <- data.frame()

for(i in 1:length(IDs)){
  x <- samps[which(samps$Individual.ID == IDs[i]), cols_same]
  dat <- rbind(dat, x[1,])
}       

write.csv(dat, "FACS_SFARI_fetalSamples_IndividualInfo.csv")


# make info columns for data that varied between samples

x <- 1:2
y <- 1:10
n <- max(length(x), length(y))
length(x) <- n                      
length(y) <- n



cols_samp_names <- c("Individual.ID", paste0("Sample_Name_", 1:7))

dat_diff <- data.frame()
for(i in 1:length(IDs)){
  x <- samps[which(samps$Individual.ID == IDs[i]),]
  Sample_Names <- t(x$Sample.Name)
  Sample_Names <- c(x$Individual.ID[1], Sample_Names)
  names(Sample_Names) <- c("Individual.ID", paste0("Sample_Name_", 1:nrow(x)))
  dat_diff <- as.data.frame(bind_rows(dat_diff, Sample_Names))
  
}

full_dat <- left_join(dat, dat_diff)

write.csv(full_dat, "FACS_SFARI_fetalSamples_IndividualandSampleInfo.csv")

# cell fractions
 
# change to standarsized codes
# note top/bottom/middle need to be changed to total, and note made in notes col for these samples

table(samps$Cell.Fraction)

samps$Cell.Fraction <- gsub("BotFraction|Middle|TopFraction|Total Nuclei", "Total", samps$Cell.Fraction)
samps$Cell.Fraction <- gsub("NeuN", "Neun_pos", samps$Cell.Fraction)
samps$Cell.Fraction <- gsub("DoubleNegative", "Double_neg", samps$Cell.Fraction)
samps$Cell.Fraction <- gsub("SATB2-", "Satb2_neg", samps$Cell.Fraction, fixed = T)
samps$Cell.Fraction <- gsub("SATB2+", "Satb2_pos", samps$Cell.Fraction, fixed = F)
samps$Cell.Fraction <- gsub("Sox10", "Sox10_pos", samps$Cell.Fraction)

cols_cell_fracs <- names(table(samps$Cell.Fraction))
dat_fracs <- as.data.frame(IDs)
dat_fracs[cols_cell_fracs] <- NA
colnames(dat_fracs)[1] <- colnames(full_dat)[1]

for(i in 1:length(IDs)){
  x <- samps[which(samps$Individual.ID == IDs[i]),]
  dat_fracs[i, "Double_neg"] <- length(which(x$Cell.Fraction == "Double_neg"))
  dat_fracs[i, "Neun_pos"] <- length(which(x$Cell.Fraction == "Neun_pos"))
  dat_fracs[i, "Satb2_neg"] <- length(which(x$Cell.Fraction == "Satb2_neg"))
  dat_fracs[i, "Satb2_pos+"] <- length(which(x$Cell.Fraction == "Satb2_pos+"))
  dat_fracs[i, "Sox10_pos"] <- length(which(x$Cell.Fraction == "Sox10_pos"))
  dat_fracs[i, "Total"] <- length(which(x$Cell.Fraction == "Total"))
}

full_dat <- left_join(full_dat, dat_fracs)

write.csv(full_dat, "FACS_SFARI_fetalSamples_final.csv")
