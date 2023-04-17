library(Seurat)
library(here)
library(sctransform)
library(glmGamPoi)
library(dplyr)


# 1. Integrate three control and three potassium deficient samples-------------------
# Load six pre-processed data set and merge
nk1 <- read_rds(here("NK1.rds"))
nk2 <- read_rds(here("NK2.rds"))
nk3 <- read_rds(here("NK3.rds"))
kd1 <- read_rds(here("KD1.rds"))
kd2 <- read_rds(here("KD2.rds"))
kd3 <- read_rds(here("KD3.rds"))

SO <- merge(nk1, y = c(nk2, nk3, kd1, kd2, kd3), add.cell.ids = c("nk1", "nk2", "nk3", "kd1", "kd2", "kd3"), project = "Diet K DCT")

# Filter low quality cells
SO <- subset(SO, subset = nFeature_RNA > 500
              & nFeature_RNA < 4000
              & nCount_RNA < 10000) 

# Normalize data with sctransform
SO <- SCTransform(SO, method = "glmGamPoi", vars.to.regress = c("nCount_RNA"), verbose = TRUE)

# Principal Component Analysis
SO <- RunPCA(SO, verbose = TRUE)
ElbowPlot(SO, ndims = 50)

# Integrate data with FindIntegrationAnchors

SO.list <- SplitObject(SO, split.by = "Rep")

SO.list <- lapply(X = SO.list, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = SO.list, nfeatures = 3000)

SO.list <- PrepSCTIntegration(object.list = SO.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = SO.list,
                                  reference = c(1, 2),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  anchor.features = features)

SO <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
SO <- RunPCA(SO, verbose = FALSE)

SO <- RunUMAP(SO, dims = 1:50)

SO <- FindNeighbors(SO, reduction = "pca", dims = 1:50)

SO <- FindClusters(SO, resolution = 0.2)
DefaultAssay(SO) <- "SCT"

# Check batch effects
DimPlot(SO,
        reduction = "umap",
        label = TRUE,
        pt.size = 0,
        label.size = 6) +
  NoLegend()

SO$Rep <- factor(x = SO$Rep, levels = c("NK1", "NK2", "NK3", "KD1", "KD2", "KD3")) 
DimPlot(SO, reduction = "umap", split.by = "Rep", ncol = 3)

# 2. Annotate clusters ------------------------------------------------

# Quality control
VlnPlot(SO, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        fill.by = "ident",
        stack = TRUE, 
        flip = TRUE,
        pt.size = 0) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, hjust = .5),
        axis.title.x = element_blank()) +
  stat_summary(fun = median,
               geom = "crossbar",
               width = 0.3,
               size = 0.2,
               position = position_dodge(width = 0.5))

FeaturePlot(SO,
            features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
            cols = c("lightgrey", "royal blue"),
            ncol = 3)

# Verification of Cluster Identity using canonical markers
canonical_markers <- c("Slc12a3",      # DCT
                       "Slc8a1",       # DCT2, CNT
                       "Slc5a12",      # PT-S1
                       "Egf",          # TAL, DCT1
                       "Umod",         # TAL, DCT1
                       "Slc12a1",      # TAL
                       "Aqp2",         # PC
                       "Slc26a4",      # IC-B
                       "Lrp2",         # PT
                       "Flt1",         # Endo
                       "Slc4a1",       # IC-A
                       "Pdgfrb",       # Peri
                       "Top2a",        # Proliferation
                       "Nphs1",        # Podo
                       "Pvalb",        # DCT1
                       "Ptprc")        # Immune

VlnPlot(SO,
        features = canonical_markers,
        stack = TRUE,
        flip = TRUE,
        pt.size = 0,
        fill.by = "ident") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = .5),
        axis.title.x = element_blank()) + 
  stat_summary(fun = median,
               geom = "crossbar",
               width = 0.3,
               size = 0.1,
               position = position_dodge(width = 0.5))

# Top five genes in each cluster to identify them
SO <- PrepSCTFindMarkers(SO, verbose = TRUE)
SO.markers <- FindSOMarkers(SO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SO.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- SO.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

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
  coord_flip()

# Name each cluster based on their identities ---- add to meta.data
SO@meta.data <- SO@meta.data %>% mutate(class = dplyr::case_when(
  seurat_clusters == 0 ~ "DCT1",
  seurat_clusters == 1 ~ "DCT1",
  seurat_clusters == 2 ~ "DCT2",
  seurat_clusters == 3 ~ "PT-S1",
  seurat_clusters == 4 ~ "TAL",
  seurat_clusters == 5 ~ "PC",
  seurat_clusters == 6 ~ "IC-B",
  seurat_clusters == 7 ~ "PT-S2",
  seurat_clusters == 8 ~ "EC",
  seurat_clusters == 9 ~ "IC-A",
  seurat_clusters == 10 ~ "Peri",
  seurat_clusters == 11 ~ "Prolif",
  seurat_clusters == 12 ~ "X",
  seurat_clusters == 13 ~ "Podo"))

# Visualize the clusters
Idents(SO) <- "class"

DimPlot(SO,
        reduction = "umap",
        label = TRUE,
        label.size = 5,
        pt.size = 0.5,
        cols = mycols) +
  NoLegend()

# 3. Subset DCT nuclei ----------------------------------------------------------

# Subset DCT1, DCT2, Proliferating population
SO.DCT <- subset(SO, idents = c("DCT1", "DCT2", "Prolif"))

# Re-analyze DCT nuclei in SO six samples
DefaultAssay(SO.DCT) <- "integrated"
SO.DCT <- RunPCA(SO.DCT, verbose = FALSE)
SO.DCT <- RunUMAP(SO.DCT, dims = 1:50)
SO.DCT <- FindNeighbors(SO.DCT, reduction = "pca", dims = 1:50)
SO.DCT <- FindClusters(SO.DCT, resolution = 0.1)
DefaultAssay(SO.DCT) <- "SCT"

# Top five genes in each cluster to identify them
SO.DCT <- PrepSCTFindMarkers(SO.DCT, verbose = TRUE)
SO.markers.DCT <- FindSOMarkers(SO.DCT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SO.markers.DCT %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- SO.markers.DCT %>% group_by(cluster) %>% top_n(5, avg_log2FC)

DotPlot(SO.DCT,
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
  RotatedAxis()

# Name each cluster based on their identities ---- add to meta.data

## Adding Annotation Metadata to DCT Clusters
SO.DCT@meta.data <- SO.DCT@meta.data %>% 
  mutate(subclass.l1 = dplyr::case_when(
    seurat_clusters == 0 ~ "DCT1",
    seurat_clusters == 1 ~ "DCT2",
    seurat_clusters == 2 ~ "Prolif"))

## Adding Annotation Metadata to DCT Clusters (all cells as DCT)
SO.DCT@meta.data <- SO.DCT@meta.data %>% 
  mutate(class = dplyr::case_when(
    seurat_clusters == 0 ~ "DCT",
    seurat_clusters == 1 ~ "DCT",
    seurat_clusters == 2 ~ "DCT"))

## This creates another metadata, 'analysis' with cell type and treatment, which will be used for later DEG analysis.
SO.DCT$analysis <- paste(SO.DCT$subclass.l1, SO.DCT$Diet, sep = "_")

# 4. Save the final seurat object
Idents(SO.DCT) <- "subclass.l1"
my_levels <- c("DCT1", "DCT2", "Prolif")
Idents(SO.DCT) <- factor(x = Idents(SO.DCT), levels = my_levels)

saveRDS(SO.DCT, "ALL_DCT.rds")

























