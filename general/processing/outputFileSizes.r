# Emma Walker
# 14/01/2022
# E.M.Walker@exeter.ac.uk

# output file sizes as csv file

setwd("/gpfs/mrc0/projects/Research_Project-MRC190311/ATACSeq/RawData/Simons/ftp1.sequencing.exeter.ac.uk/V0161/01_raw_reads/")

dat <- as.data.frame(file.info(list.files(".")))
dat <- dat %>% dplyr::select(size)
dat$size2 <- apply(dat, 1, utils:::format.object_size, "auto")

write.csv(dat, "../../../RawDataFileSizes.csv")




