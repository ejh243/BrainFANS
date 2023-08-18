##---------------------------------------------------------------------#
##
## Title: Compute Correlation Pairwise
##
## Purpose of script: Compute corr sites in common pairwise
##
## Author: Jessica Shields
##
## Date Created: 2022-06-08
##
##---------------------------------------------------------------------#
tissueFile<-list.files(paste0(methylDir, '/QCOutput/'), pattern = '.chr1.bg') 

for (i in 2:length(tissueFile)){
  covMethyl<-fread(paste(methylDir, 'QCOutput', tissueFile[i], sep = '/'), sep = '\t') %>%
    .[, c(2, 4, 8)] %>%
    pivot_wider(names_from = V4, values_from = V8)
  
  n<-length(colnames(covMethyl))-1
  
  # create empty matrix
  corrMat<-NULL
  corrMat<- matrix(data = NA, nrow = n, ncol = n-1)
  rownames(corrMat)<-colnames(covMethyl)[2:length(colnames(covMethyl))]
  colnames(corrMat)<-colnames(covMethyl)[2:(length(colnames(covMethyl))-1)]
  
  # populate matrix
  for (y in 1:(n-1)) {
    for (x in (y+1):n) {
      # fill in correlation matrix with site percentages from the imported data  
      corrMat[x,y] <- cor(covMethyl[(x+1)], covMethyl[(y+1)],  method = "pearson", use = "complete.obs")
    }
  }
  tissue<-gsub(".chr1.bg", "", tissueFile[i])
  print(paste0(methylDir, '/QCOutput/',tissue,'.corr.genCov.qc'))
  # write to file
  #write.table(corrMat, paste0(methylDir, '/QCOutput/',tissue,'.corr.genCov.qc'), sep='\t', quote = FALSE)
}

i=1
covMethyl<-fread(paste(methylDir, 'QCOutput', tissueFile[i], sep = '/'), sep = '\t') %>%
  .[, c(2, 4, 8)] %>%
  pivot_wider(names_from = V4, values_from = V8)

n<-length(colnames(covMethyl))-1

covMethyl<-covMethyl[,order(colnames(covMethyl[1:ncol(covMethyl)]))]

# create empty matrix
corrMat<-NULL
corrMat<- matrix(data = NA, nrow = n, ncol = n-1)
rownames(corrMat)<-colnames(covMethyl)[1:length(colnames(covMethyl))-1]
colnames(corrMat)<-colnames(covMethyl)[1:(length(colnames(covMethyl))-2)]

# populate matrix
for (y in 1:(n-1)) {
  for (x in (y+1):n) {
    # fill in correlation matrix with site percentages from the imported data  
    corrMat[colnames(covMethyl[x]),
            colnames(covMethyl[y])] <- cor(covMethyl[(x)],covMethyl[(y)],  method = "pearson", use = "complete.obs")
  }
}
tissue<-gsub(".chr1.bg", "", tissueFile[i])
print(paste0(methylDir, '/QCOutput/',tissue,'.corr.genCov.qc'))
# write to file
#write.table(corrMat, paste0(methylDir, '/QCOutput/',tissue,'.corr.genCov.qc'), sep='\t', quote = FALSE)

#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#
pheno<-read.table(sampleSheet, header = TRUE, sep = ',', stringsAsFactors = FALSE)
fraction<-unique(pheno$fraction)[1]

#filter by region
pheno<-pheno[which(gsub('[()]', '', pheno$tissue) %in% gsub('_', ' ', tissue)),]

## SAME cell type
#filter by fraction 
samples<-pheno[which(pheno$fraction %in% fraction),1] 

n<-length(samples)
len<-(n-1)*n/2

pdf(paste0(qcDir, "/pairwiseCorr_",tissue,"_",fraction,".pdf"), width = 10, height = 5, onefile = TRUE)
for (y in 1:(n-1)) {
  for (x in (y+1):n) {
    print(paste0(samples[y],' compared to ', samples[x]))
    print(covMethyl[,c('V2', samples[x], samples[y])] %>% 
      drop_na() %>%
      sample_n(., 15000) %>%
      ggplot(., aes_string(x= names(.)[2], y = names(.)[3]))+
      geom_point()+
 #     labs(title = 'Brain region BA24 - sites in common pairwise\nsame individ, different cell type')+
      theme_bw()+
      geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
      annotate('text', x = 10, y = 96, label = paste0('r = ',round(corrMat[samples[x],samples[y]], 2))))
  }
}
dev.off()


## SAME individual, DIFFERENT cell type
samples<-pheno[which(pheno$fraction %in% fraction),1] 

#pull out duplicated individualIDs
pheno[pheno$individualID %in% names(which(table(pheno$individualID) > 1)),]


ids<-names(which(table(pheno$individualID) > 1))
for (x in 1:length(ids)){
  samples<-filter(pheno, individualID == ids[x])[,1]
  n<-length(samples)
  len<-(n-1)*n/2
  print(ids[x])
  pdf(paste0(qcDir, "/pairwiseCorr_",tissue,"_",ids[x],".pdf"), width = 10, height = 5, onefile = TRUE)
  for (y in 1:(n-1)) {
    for (x in (y+1):n) {
      print(paste0(samples[y],' compared to ', samples[x]))
      print(covMethyl[,c('V2', samples[x], samples[y])] %>% 
              drop_na() %>%
              sample_n(., 15000) %>%
              ggplot(., aes_string(x= names(.)[2], y = names(.)[3]))+
              geom_point()+
              #     labs(title = 'Brain region BA24 - sites in common pairwise\nsame individ, different cell type')+
              theme_bw()+
              geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
              annotate('text', x = 10, y = 96, label = paste0('r = ',round(corrMat[samples[x],samples[y]], 2))))
    }
  }
  dev.off()
}

## DIFFERENT individual, DIFFERENT cell type






pdf(paste0(qcDir, "/pairwiseCorr_",fraction,".pdf"), width = 10, height = 10)
samples=c("SRR5343780", "SRR5343781")
covMethyl[,c('V2', samples[1], samples[2])] %>% 
  drop_na() %>%
  sample_n(., 15000) %>%
  ggplot(., aes(x=SRR5343780, y =SRR5343781))+
  geom_point()+
  labs(title = 'Brain region BA24 - sites in common pairwise\nsame individ, different cell type', 
       x = 'SRR5343780 - 5248_BA24_neg', y = 'SRR5343781 - 5248_BA24_pos')+
  theme_bw()+
  geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
  annotate('text', x = 10, y = 96, label = paste0('r = ',corrMat[samples[2],samples[1]]))


samples=c("SRR5343780", "SRR5343788")
covMethyl[,c('V2', samples[1], samples[2])] %>% 
  drop_na() %>%
  sample_n(., 15000) %>%
  ggplot(., aes_string(x= names(.)[2], y = names(.)[3]))+
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
  labs(title = 'Brain region BA24 - sites in common pairwise\ndifferent individ, different cell type', 
       x = 'SRR5343780 - 5248_BA24_neg', y = 'SRR5343788 - 5284_BA24_pos')+
  theme_bw()+
  annotate('text', x = 10, y = 96, label = paste0('r = ',corrMat[samples[2],samples[1]]))

samples=c("SRR5343780", "SRR5343795")
covMethyl[,c('V2', samples[1], samples[2])] %>% 
  drop_na() %>%
  sample_n(., 15000) %>% data.frame(.) %>%
  ggplot(., aes_string(x= names(.)[2], y = names(.)[3]))+
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
  labs(title = 'Brain region BA24 - sites in common pairwise\ndifferent individ, same cell type', 
       x = 'SRR5343780 - 5248_BA24_neg', y = 'SRR5343795 - 5552_BA24_neg')+
  theme_bw()+
  annotate('text', x = 15, y = 96, label = paste0('r = ',corrMat[samples[2],samples[1]]))


samples=c("SRR5343780", "SRR5343802")
covMethyl[,c('V2', samples[1], samples[2])] %>% 
  drop_na() %>%
  sample_n(., 15000) %>% data.frame(.) %>%
  ggplot(., aes_string(x= names(.)[2], y = names(.)[3]))+
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
  labs(title = 'Brain region BA24 - sites in common pairwise\ndifferent individ, same cell type', 
       x = 'SRR5343780 - 5248_BA24_neg', y = 'SRR5343795 - 5552_BA24_neg')+
  theme_bw()+
  annotate('text', x = 15, y = 96, label = paste0('r = ',corrMat[samples[2],samples[1]]))

samples=c("SRR5343795", "SRR5343802")
covMethyl[,c('V2', samples[1], samples[2])] %>% 
  drop_na() %>%
  sample_n(., 15000) %>% data.frame(.) %>%
  ggplot(., aes_string(x= names(.)[2], y = names(.)[3]))+
  geom_point()+
  geom_smooth(method = 'lm', se = FALSE, colour = 'black')+
  labs(title = 'Brain region BA24 - sites in common pairwise\ndifferent individ, same cell type', 
       x = 'SRR5343780 - 5248_BA24_neg', y = 'SRR5343795 - 5552_BA24_neg')+
  theme_bw()+
  annotate('text', x = 15, y = 96, label = paste0('r = ',round(corrMat[samples[2],samples[1]], 2)))

pheno$fraction %in% tissue        
