
setwd("C:/Users/pgand/OneDrive/Documents/Bioinformatics Journey/")


# Install Packages
#install.packages("Seurat")
#install.packages("tidyverse")

# Install Libraries
library('Seurat')
library('tidyverse')

# Install Raw Data
hum_glioblastoma_obj <- Read10X_h5(filename = 'Human Glioblastoma Multiforme Transcriptome Analysis/Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.h5')

# Initializing the Seurat object with raw (non-normalized data)
hum_glioblastoma_obj <- CreateSeuratObject(counts = hum_glioblastoma_obj, project="brain", min.cell=3, min.features=200)
str(hum_glioblastoma_obj)
hum_glioblastoma_obj

# Now that the the data is loaded, we have our Seurat Object, and know number the number of features (24381) and samples (5960).
# We can start the steps of analyzing the data.

# 1. QC and selecting cells for further analysis
hum_glioblastoma_obj[['percent.mt']] <- PercentageFeatureSet(hum_glioblastoma_obj, pattern = "^MT-")
View(hum_glioblastoma_obj@meta.data)
head(hum_glioblastoma_obj, 10)

# 1b. Visualization of QC metrics
VlnPlot(hum_glioblastoma_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 1c. FeatureScatter: Feature - Feature relationship
plot1 <- FeatureScatter(hum_glioblastoma_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hum_glioblastoma_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
plot1 + plot2

# 2. Filter out low quality cells
hum_glioblastoma_obj <- subset(hum_glioblastoma_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hum_glioblastoma_obj

# 3. Normalization
hum_glioblastoma_obj <- NormalizeData(hum_glioblastoma_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Identify highly variable features (Cell-to-Cell Variation) 
hum_glioblastoma_obj <- FindVariableFeatures(hum_glioblastoma_obj, selection.method = "vst", nfeatures = 2000)

# 4b. Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hum_glioblastoma_obj), 10)
top10

# 4c. plot variable features with and without labels
plot1 <- VariableFeaturePlot(hum_glioblastoma_obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling the data (linear transformation)
all.genes <- rownames(hum_glioblastoma_obj)
hum_glioblastoma_obj <- ScaleData(hum_glioblastoma_obj, features = all.genes)

# 6. Perform linear dimensional reduction
hum_glioblastoma_obj <- RunPCA(hum_glioblastoma_obj, features = VariableFeatures(object = hum_glioblastoma_obj))

# 6b. Examine and visualize PCA results a few different ways
print(hum_glioblastoma_obj[["pca"]], dims = 1:5, nfeatures = 5)

# 6c.
DimHeatmap(hum_glioblastoma_obj, dims = 1, cells = 500, balanced = TRUE)

# 6d.
VizDimLoadings(hum_glioblastoma_obj, dims = 1:2, reduction = "pca")

# 6e.
DimPlot(hum_glioblastoma_obj, reduction = "pca")

# 7. Determine the dimensionality of the dataset
ElbowPlot(hum_glioblastoma_obj)

# 7b.
hum_glioblastoma_obj <- JackStraw(hum_glioblastoma_obj, num.replicate = 100)

# 7c.
hum_glioblastoma_obj <- ScoreJackStraw(hum_glioblastoma_obj, dims = 1:20)

# 7d.
JackStrawPlot(hum_glioblastoma_obj, dims = 1:15)

# 8. Cluster the cells
hum_glioblastoma_obj <- FindNeighbors(hum_glioblastoma_obj, dims = 1:15)

# 8b.
hum_glioblastoma_obj <- FindClusters(hum_glioblastoma_obj, resolution = c(0.1, 0.3, 0.5, .8, 1))
View(hum_glioblastoma_obj@meta.data)

DimPlot(hum_glioblastoma_obj, group.by = "RNA_snn_res.1", label = TRUE)

# 9. Run non-linear dimensional reduction (UMAP/tSNE)
hum_glioblastoma_obj <- RunUMAP(hum_glioblastoma_obj, dims = 1:15)

# 9b. individual clusters
DimPlot(hum_glioblastoma_obj, reduction = "umap")

# 10. Finding differently expressed features (cluster biomarkers)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(hum_glioblastoma_obj, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# 10b. find all markers distinguishing cluster 7 from clusters 8 and 9
cluster7.markers <- FindMarkers(hum_glioblastoma_obj, ident.1 = 7, ident.2 = c(8, 9), min.pct = 0.25)
head(cluster7.markers, n = 7)

# 10c. find markers for every cluster compared to all remaining cells, report only the positive
hum_glioblastoma_obj.markers <- FindAllMarkers(hum_glioblastoma_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hum_glioblastoma_obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


# 10d
FeaturePlot(hum_glioblastoma_obj, features = c("GFAP", "SERPINA3", "QDPR", "C11orf96", "GPR183", "HLA-DRB1", "FCGBP", "C1QC",
                              "TIMP1", "S100A9"))

# 10e
cluster0.markers <- FindMarkers(hum_glioblastoma_obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# 10f
hum_glioblastoma_obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(hum_glioblastoma_obj, features = top10$gene) + NoLegend()

