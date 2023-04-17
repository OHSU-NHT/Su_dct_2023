
library(Seurat)
library(here)
library(dplyr)

# 1. Sub population of DCT1 -------------------------------------

# Subset DCT1 cells
SO_dct <- readRDS(here("CONTROL_DCT.rds"))
Idents(SO_dct) <- "subclass.l1"

SO <- subset(SO_dct, idents = c("DCT1"))
DefaultAssay(SO) <- "integrated"

# Re-analyze DCT1 cells
SO <- RunPCA(SO, verbose = FALSE)
ElbowPlot(SO, ndims = 50)
SO <- RunUMAP(SO, dims = 1:50)
SO <- FindNeighbors(SO, reduction = "pca", dims = 1:50)
SO <- FindClusters(SO, resolution = 0.2)
DefaultAssay(SO) <- "SCT"

# Find top 5 genes in each cluster
SO <- PrepSCTFindMarkers(SO, verbose = TRUE)
ALL.markers.DCT1.sub <- FindAllMarkers(SO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ALL.markers.DCT1.sub %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- ALL.markers.DCT1.sub %>% group_by(cluster) %>% top_n(5, avg_log2FC)

DotPlot(SO,
        features = unique(top5$gene),
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        scale.max = 100,
        scale.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  RotatedAxis() +
  scale_y_discrete(limits = c("3", "2", "1", "0"))

#Adding Annotation Metadata to DCT Clusters
SO@meta.data <- SO@meta.data %>% 
  mutate(subclass.l2 = dplyr::case_when(
    seurat_clusters == 0 ~ "DCT1-A",
    seurat_clusters == 1 ~ "DCT1-D",
    seurat_clusters == 2 ~ "DCT1-B",
    seurat_clusters == 3 ~ "DCT1-C"))

Idents(SO) <- "subclass.l2"
my_levels <- c("DCT1-A", "DCT1-B", "DCT1-C", "DCT1-D")
Idents(SO) <- factor(x = Idents(SO), levels = my_levels)

saveRDS(SO, file = "CONTROL_DCT1.rds")

# 2. Sub population of DCT2 -------------------------------------

# Subset DCT2 cells
SO_dct <- readRDS(here("CONTROL_DCT.rds"))
Idents(SO_dct) <- "subclass.l1"

SO <- subset(SO_dct, idents = c("DCT2"))
DefaultAssay(SO) <- "integrated"

# Re-analyze DCT2 cells
SO <- RunPCA(SO, verbose = FALSE)
ElbowPlot(SO, ndims = 50)
SO <- RunUMAP(SO, dims = 1:10)
SO <- FindNeighbors(SO, reduction = "pca", dims = 1:10)
SO <- FindClusters(SO, resolution = 0.15)
DefaultAssay(SO) <- "SCT"

# Find top 5 genes in each cluster
SO <- PrepSCTFindMarkers(SO, verbose = TRUE)
ALL.markers.DCT1.sub <- FindAllMarkers(SO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ALL.markers.DCT1.sub %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- ALL.markers.DCT1.sub %>% group_by(cluster) %>% top_n(5, avg_log2FC)

DotPlot(SO,
        features = unique(top5$gene),
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        scale.max = 100,
        scale.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  RotatedAxis() +
  scale_y_discrete(limits = c("1", "0"))

#Adding Annotation Metadata to DCT Clusters
SO@meta.data <- SO@meta.data %>% 
  mutate(subclass.l2 = dplyr::case_when(
    seurat_clusters == 0 ~ "DCT2-beta",
    seurat_clusters == 1 ~ "DCT2-alpha"))

Idents(SO) <- "subclass.l2"
my_levels <- c("DCT2-alpha", "DCT2-beta")
Idents(SO) <- factor(x = Idents(SO), levels = my_levels)

saveRDS(SO, file = "CONTROL_DCT2.rds")
