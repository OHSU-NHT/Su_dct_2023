
library(Seurat)
library(here)
library(DESeq2)
library(dplyr)
library(openxlsx)

# 1. Comparing DCT1 and DCT2 in control samples ----------------------------------

# Load CONTROL_DCT seurat object
SO <- readRDS(here("CONTROL_DCT.rds"))
Idents(SO) <- "subclass.l1"

# Calculate DEGs comparing DCT1 and DCT2 using test DESeq2
DCT1vsDCT2 <- FindMarkers(SO, ident.1 = "DCT1", ident.2 = "DCT2", test.use = "DESeq2", verbose = TRUE)

# 2. Computing Proliferating cell DEGs in control samples ------------------------ 
Prolif <- FindMarkers(SO, ident.1 = "Prolif", ident.2 = c("DCT1", "DCT2"), test.use = "DESeq2", verbose = TRUE)

# 3. DEGs of DCT1 and DCT2 subpopulation in control samples ----------------------------------

# Load CONTROL_DCT1 seurat object
SO <- readRDS(here("CONTROL_DCT1.rds"))
Idents(SO) <- "subclass.l2"

DCT1a <- FindMarkers(SO,
                        ident.1 = "DCT1-A",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

DCT1b <- FindMarkers(SO,
                        ident.1 = "DCT1-B",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

DCT1c <- FindMarkers(SO,
                        ident.1 = "DCT1-C",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

DCT1d <- FindMarkers(SO,
                        ident.1 = "DCT1-D",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

# Load CONTROL_DCT2 seurat object
SO <- readRDS(here("CONTROL_DCT2.rds"))
Idents(SO) <- "subclass.l2"

DCT2x.vs.y <- FindMarkers(SO,
                        ident.1 = "DCT2-X",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

# 4. Save all DEG lists -----------------------------------------------------
# data.frame
save(DCT1vsDCT2, Prolif, DCT1a, DCT1b, DCT1c, DCT1d, DCT2x.vs.y, file = here("DCT DEGs.RData"))
save(DCT1vsDCT2, Prolif, DCT1a, DCT1b, DCT1c, DCT1d, DCT2x.vs.y, file = here("DCT DEGs.rda"))

# excel file
dataset_names <- list('DCT1vsDCT2' = DCT1vsDCT2, 'Prolif' = Prolif, 'DCT1a' = DCT1a, 'DCT1b' = DCT1b, 'DCT1c' = DCT1c, 'DCT1d' = DCT1d, 'DCT2x.vs.y' = DCT2x.vs.y)
write.xlsx(dataset_names, file = 'DCT DEGs.xlsx', rowNames=TRUE)
