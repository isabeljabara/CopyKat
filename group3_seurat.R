library(dplyr)
library(Seurat)
library(patchwork)
library(Cairo)

#load patients data, returns UMI count matrix
tnbc.data <- Read10X(data.dir = "/Users/kristenfrombach/Desktop/N1B_RAW/")
#make a seurat object with raw data, though i think this has already been
#normalized by authors
tnbc <- CreateSeuratObject(counts = tnbc.data, project = "tnbcpt1", min.cells = 3, min.features = 200)
#tnbc = An object of class Seurat 
#19141 features (genes) across 3634 samples (individual cells) within 1 assay 
#Active assay: RNA (19141 features, 0 variable features)

#QC - REMOVING UNWANTED CELLS
#labeling mitochondrial
tnbc[["percent.mt"]] <- PercentageFeatureSet(tnbc, pattern = "^MT-")
#visualizing # of features, RNA count, mito %
#pt.size = 0 because data too dense, can't see distribution (can take out if want full)
VlnPlot(tnbc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
#WEIRD ERROR FOR VLNPLOT; WILL NEED TO MAYBE UPDATE AND QUIT R???
#plot1 <- FeatureScatter(tnbc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(tnbc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#given plots from above, mito cutoff 20%, feature min 200,max 25000 
tnbc <- subset(tnbc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10)
#removed unwanted cell

#NORMALIZE
tnbc <- NormalizeData(tnbc, normalization.method = "LogNormalize", scale.factor = 10000)

#ID VARIABLE FEATURES
tnbc <- FindVariableFeatures(tnbc, selection.method = "vst", nfeatures = 2000)
#ID 20 most highly variable genes -- hoping our proteins of interest there?
#top20 <-head(VariableFeatures(tnbc), 20)
#plot these variable features with and w/o labels
#plot1 <- VariableFeaturePlot(tnbc)
#plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#plot1 + plot2

#SCALING DATA
#shifts expression of each gene so mean across cell is 0
#scales so that variance of each gene across cells is 1
#results stored in tnbc[["RNA']]@scale.data
all.genes <- rownames(tnbc)
tnbc <- ScaleData(tnbc, features = all.genes)

#LINEAR DIMENSIONAL REDUCTION - PCA
tnbc <- RunPCA(tnbc, features = VariableFeatures(object = tnbc))
#visualize PCA results
#print(tnbc[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(tnbc, dims = 1:2, reduction = "pca")
#DimPlot(tnbc, reduction = "pca")
#DimHeatmap(tnbc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(tnbc, dims = 1:15, cells = 500, balanced = TRUE)
#determine dimensionality
ElbowPlot(tnbc)
#could choose 15-20,will go with 15 for now but can also try 20 later

#CLUSTERING
#resolution sets 'granularity', .4-1.2 typically good for 3k cell sample
#(ours ~3.5k so close), optimal resolution increase for larger datasets
tnbc <- FindNeighbors(tnbc, dims = 1:15)
tnbc <- FindClusters(tnbc, resolution = 0.4)
#Found 12 communities

#NON LINEAR DIM REDUCTION - UMAP/tSNE
tnbc <- RunUMAP(tnbc, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(tnbc, reduction = "umap")

#FIND CLUSTER BIOMARKERS - hopefully can ID clusters
# find all markers of cluster 2
# cluster2.markers <- FindMarkers(tnbc, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(tnbc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
tnbc.markers <- FindAllMarkers(tnbc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster7.markers <- FindMarkers(tnbc, ident.1 = 6, min.pct = 0.25)
head(cluster7.markers, n = 10)
cluster0.markers <- FindMarkers(tnbc, ident.1 = 8, min.pct = 0.25)
head(cluster0.markers, n = 10)
tnbc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
tnbc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(tnbc, features = top10$gene) + NoLegend()
FeaturePlot(tnbc, features = c("EPCAM", "ITGA6", "CD3D", "PTPRC", "FCGR3A", "CSF1R"))
FeaturePlot(tnbc, features = c("TACSTD2", "SLC39A6", "ERBB2", "ERBB3", "PTK7", "ROR2"))
new.cluster.ids <- c("Epithelial 1", "Undifferentiated", "Keratinocytes", "Fibroblasts", "Epithelial 2", "Immune")
names(new.cluster.ids) <- levels(tnbc)
tnbc <- RenameIdents(tnbc, new.cluster.ids)
DimPlot(tnbc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# TACSTD2 = TROP2
# SLC39A6 = LIV1
# ERBB2 = HER2
# ERBB3 = HER3
# PTK7 = PTK7
# ROR2 = ROR2
