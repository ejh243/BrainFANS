#https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
#set working directory and load in R packages 
path = 
setwd(path)
library(Seurat)
library(dplyr)
library(patchwork)

# Load the dataset
datapath = 
scdata.data <- Read10X(data.dir = datapath)

# Initialize the Seurat object with the raw (non-normalized) data. Minimum of 200 genes per cell
scdata <- CreateSeuratObject(counts = scdata.data, project = "X", min.cells = 3, min.features = 200)
scdata

#Look at the size of the matrix

dense.size <- object.size(as.matrix(scdata.data))
dense.size

sparse.size <- object.size(scdata.data)
sparse.size

dense.size/sparse.size

# Use the [[ operator to add columns to object metadata. Add % mitochondrial genes to the QC stats
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(scdata@meta.data, 5)

# Save QC metrics into a table
write.table(scdata@meta.data, "scdata_QC_metrics.txt", sep="\t")


#Visualize QC metrics as a violin plot and save as pdf
#make a pdf, in this location, with this name
#create a violin plot, which will auomatically be added to the made pdf.options
#dev.off means close the pdf so that nothing else can be added to it 
pdf("scdata_QCmetric.pdf")
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#Use FeatureScatter to visualize feature-feature relationships- particularly important for QC is the feature counts vs Mitochondrial genes 

pdf("scdata_ncount_mt.pdf")
plot1 <- FeatureScatter(scdata, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = "darkblue")
plot1
dev.off()

pdf("scdata_ncount_featurecount.pdf")
plot2 <- FeatureScatter(scdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5, cols = "darkblue")
plot2
dev.off()

pdf("scdata_featurecount_mt.pdf")
plot3 <- FeatureScatter(scdata, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = "darkblue")
plot3
dev.off()

#Filter out cells that have unique feature counts over 2,500 or less than 200 and that have >5% mitochondrial counts

scdata <- subset(scdata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
scdata

#Normalise the data: normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.

scdata <- NormalizeData(scdata, normalization.method = "LogNormalize", scale.factor = 10000)

#Calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). Return 2000 features per dataset

scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scdata), 10)
top10

# plot variable features with and without labels
pdf("scdata_Variablefeature.pdf")
plot1 <- VariableFeaturePlot(scdata)
plot1
dev.off()

pdf("scdata_Variable_Feature_Labelled.pdf")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge=0, ynudge=0)
plot2
dev.off()

#Apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA

all.genes <- rownames(scdata)
scdata <- ScaleData(scdata, features = all.genes)

#You can perform scaling only on the previously identified variable features (2,000 by default). To do this, omit the features argument

#Can 'regress out' variation caused by things like cell cycle and MT contamination- advice to use new feature called sctransform
#MAY WANT TO MISS OUT THIS STEP
scdata_clean <- ScaleData(scdata, vars.to.regress = "percent.mt")

#Perform PCA on the scaled data- uses only varibale features identified previously

scdata_clean_pca <- RunPCA(scdata_clean, features = VariableFeatures(object = scdata_clean))

pdf("scdata_clean_PCA.pdf")
DimPlot(scdata_clean_pca, reduction = "pca")
dev.off()

#we also went back and used the uncleaned data just to compare with and without MT. so we also had to create scdata_pca(not clean)
scdata_pca <- RunPCA(scdata, features = VariableFeatures(object = scdata))

pdf("scdata_PCA.pdf")
DimPlot(scdata_pca, reduction = "pca")
dev.off()

#then back to Gina's process.
pdf("scdata_PC_heatmaps.pdf")
DimHeatmap(scdata_clean_pca, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
#it might be worth seeing if we can introduce a plot to show which PCAs contribute to the variability and then subsequently only use these PCAs

pdf("elbowplot_pca.pdf")
ElbowPlot(scdata_clean_pca, ndims = 50)
dev.off()

#Choosing which PCs to use in subsequent analysis: resampling test inspired by the JackStraw procedure. Randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of feature scores, and repeat this procedure. We identify 'significant' PCs as those who have a strong enrichment of low p-value features. 
scdata_jackstraw <- JackStraw(scdata_clean_pca, num.replicate = 100)
scdata_scorejackstraw <- ScoreJackStraw(scdata_jackstraw, dims = 1:20)

pdf("scdata_scorejackstraw.pdf")
JackStrawPlot(scdata_scorejackstraw, dims = 1:15)
dev.off()

#Alternative approach to choose PCs: Use an elbow plot to rank principle components based on the percentage of variance explained by each one 
#didn't see this, turns out Gina already suggested it
pdf("scdata_elbow_plot.pdf")
ElbowPlot(scdata)
dev.off()

#Clustering cells
#1)Construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).

scdata_neighbor <- FindNeighbors(scdata_scorejackstraw, dims = 1:15)

#2) Apply modularity optimization techniques such as the Louvain algorithm (default) or SLM to iteratively group cells together, with the goal of optimizing the standard modularity function.
# resolution parameter that sets the 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters.Setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells

scdata_cluster <- FindClusters(scdata_neighbor, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(scdata_cluster), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
#As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.

scdata_uMAP <- RunUMAP(scdata_cluster, dims = 1:15)

pdf("scdata_UMAP.pdf")
DimPlot(scdata_uMAP, reduction = "umap")
dev.off()

scdata_tSNE <- RunTSNE(scdata_cluster, dims = 1:15)

pdf("scdata_tsne.pdf")
tsne <- DimPlot(scdata_tSNE, reduction = "tsne")
tsne
dev.off()

# Find markers for every cluster compared to all remaining cells, report only the positive ones
scdata.markers <- FindAllMarkers(scdata_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scdata.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#Tests for differential expression

#ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect).
cluster1.markers <- FindMarkers(scdata_cluster, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(cluster1.markers,"cluster1.markers.csv")

#Visulaising marker exression
#FeaturePlot (visualizes feature expression on a tSNE or PCA plot)


pdf("scdata_uMAP_Marker_genes.pdf", w=12, h=12)
FeaturePlot(scdata_uMAP, features = c("SOX10", "RBFOX3", "SYT1", "GAD2", "AQP4", "PLP1", "P2RY12", "PDGFRA", "LAMA2", "PDGFRB", "PTPRC", "CUX2", "RORB"))
dev.off()


pdf("scdata_tsNE_Marker_genes.pdf", w=12, h=12)
FeaturePlot(object = scdata_tSNE, features = c("SOX10", "RBFOX3", "SYT1", "GAD2", "SLC17A7", "AQP4", "PLP1", "P2RY12", "PDGFRA", "LAMA2", "PDGFRB", "CLDN5", "PTPRC"), reduction = "tsne", ncol = 3, pt.size = 1)
dev.off()

#DoHeatmap generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pdf("scdata_Top10_heatmap.pdf")
top10 <- scdata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(scdata_cluster, features = top10$gene) + NoLegend()
dev.off()

#Assigning cell type identity to clusters (using canonical markers to match cell clusters to identity)

new.cluster.ids <- c("1", "2", "3", "4", "5", "6")
names(new.cluster.ids) <- levels(scdata_tSNE)
scdata_idents <- RenameIdents(scdata_tSNE, new.cluster.ids)

pdf("scdata_tSNE_Cluster_labelled.pdf")
DimPlot(scdata_idents, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()
dev.off()
