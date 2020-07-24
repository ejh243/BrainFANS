## Combine two datasets : Sox10 + NeunPos  (https://satijalab.org/seurat/v3.1/merge_vignette.html)

#set working directory and load in R packages
path = 
setwd(path)
library(Seurat)
library(dplyr)
library(patchwork)

# Load the datasets
sox10datapath = 
neundatapath = 
  

sox10.data <- Read10X(data.dir = sox10datapath)
sox10 <- CreateSeuratObject(counts = sox10.data, project = "Sox10", min.cells = 3, min.features = 200)
sox10

neunpos.data <- Read10X(data.dir = neundatapath)
neunpos <- CreateSeuratObject(counts = neunpos.data, project = "NeunPos", min.cells = 3, min.features = 200)
neunpos

# Combine the datasets
scdata.combined <- merge(x = sox10, y = neunpos, add.cell.ids = c("sox10", "neunpos"), project = "CombinesSortedCells")
scdata.combined

head(colnames(sox10))
head(colnames(neunpos))
head(colnames(scdata.combined))
table(scdata.combined$orig.ident)


## merge by default combines raw data. We can normalise the datasets then combine the normalised datasets. 
sox10 <- NormalizeData(sox10)
neunpos <- NormalizeData(neunpos)
scdata.normalized <- merge(x = sox10, y = neunpos, add.cell.ids = c("sox10", "neunpos"), project = "CombinesSortedCells",
                         merge.data = TRUE)
GetAssayData(scdata.combined)[1:30, 1:30]
GetAssayData(scdata.normalized)[1:30, 1:30]

