
library(Seurat)
library(DoubletFinder)
library(here)
library(SoupX)
library(dplyr)


# 1. Ambient RNA removal using SoupX ------------------------------------------------------------------------------
# Load raw and filtered aligned data and estimate soup profile
tod = Seurat::Read10X(here("NK1", "raw_feature_bc_matrix")) 
toc = Seurat::Read10X(here("NK1", "filtered_feature_bc_matrix"))
sc = SoupChannel(tod,toc)

#Make the Seurat object from the filtered control data
SO <- Read10X(here("NK1", "filtered_feature_bc_matrix"))  
SO <- CreateSeuratObject(counts = SO, project = "Diet K DCT")  

#Cluster the cells with Seurat
SO <- SCTransform(SO, verbose = F)
SO <- RunPCA(SO, verbose = F)
SO <- RunUMAP(SO, dims = 1:30, verbose = F)
SO <- FindNeighbors(SO, dims = 1:30, verbose = F)
SO <- FindClusters(SO, verbose = T)

meta_SO <- SO@meta.data
umap_SO <- SO@reductions$umap@cell.embeddings

clusters <- setNames(meta_SO$seurat_clusters, rownames(meta_SO))

#Sanity check, they should be equal
length(clusters)
nrow(sc$metaData)

sc <- setClusters(sc, clusters)
sc <- setDR(sc, umap_SO)

#Estimate rho
sc = autoEstCont(sc)
#Clean the data
SO_out = adjustCounts(sc)

#Create a new Seurat Object out of the cleaned data
SO <- CreateSeuratObject(SO_out)


# 2. Doublet removal ------------------------------------------------------------------------------

# Add metadata
SO <- AddMetaData(object = SO, metadata = "NK1", col.name = "Rep") 
SO <- AddMetaData(object = SO, metadata = "NK", col.name = "Diet") 
SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^mt-")
SO

# QC and Filtering (low quality cells)
VlnPlot(SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
SO.f <- subset(SO, subset = nFeature_RNA < 5000 
               & nCount_RNA < 16000)
VlnPlot(SO.f, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
SO.f

# Pre-process standard workflow
SO.f <- NormalizeData(object = SO.f)
SO.f <- FindVariableFeatures(object = SO.f)
SO.f <- ScaleData(object = SO.f)
SO.f <- RunPCA(object = SO.f)
ElbowPlot(SO.f, 50)

# PCs between 30-40
SO.f <- FindNeighbors(object = SO.f, dims = 1:40)
SO.f <- FindClusters(object = SO.f)
SO.f <- RunUMAP(object = SO.f, dims = 1:40)



# Calculate each combination of pN and pK
sweep.res.list_SO.f <- paramSweep_v3(SO.f, PCs = 1:50, sct = FALSE) 


#Summarize each combination of pN and pK
sweep.stats_SO.f <- summarizeSweep(sweep.res.list_SO.f, GT = FALSE) 

#Select the pK that corresponds to max bcmvn to optimize doublet detection
bcmvn_SO.f <- find.pK(sweep.stats_SO.f)
pK <- bcmvn_SO.f %>% 
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 

#See pK in the Values Environment
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate 
annotations <- SO.f@meta.data$seurat_clusters  

homotypic.prop <- modelHomotypic(annotations)           
homotypic.prop

# 10X Multiplet Rate Table https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-

nExp_poi <- round(0.08*nrow(SO.f@meta.data)) 
nExp_poi
nExp_poi_adj <- round(nExp_poi*(1-homotypic.prop))


## Doublet Finder

SO.f <- doubletFinder_v3(SO.f,
                         PCs = 1:20,
                         pN = 0.25,
                         pK = pK,
                         nExp = nExp_poi_adj,
                         reuse.pANN = FALSE, sct = FALSE)
colnames(SO.f@meta.data)[9] <- "pANN"
colnames(SO.f@meta.data)[10] <- "DF.class"
head(SO.f@meta.data)
table(SO.f@meta.data$DF.class)


## Subset singlets

SO.f_singlets <- subset(SO.f, DF.class == "Singlet")


# 3. Remove mitochondrial genes ------------------------------------------------------------------------------

SO.f_singlets <- SO.f_singlets[!grepl("^mt-", rownames(SO.f_singlets)), ]

# Save the final dataset for integration
saveRDS(SO.f_singlets, here("NK1.rds")) 

## Process all other samples following the same steps above.


