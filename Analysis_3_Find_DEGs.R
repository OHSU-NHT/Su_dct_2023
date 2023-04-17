
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
DCT1vsDCT2_NK <- FindMarkers(SO, ident.1 = "DCT1", ident.2 = "DCT2", test.use = "DESeq2", verbose = TRUE)

# 2. Computing Proliferating cell DEGs in control samples ------------------------ 
Prolif_NK <- FindMarkers(SO, ident.1 = "Prolif", ident.2 = c("DCT1", "DCT2"), test.use = "DESeq2", verbose = TRUE)

# 3. Computing potassium-deficiency-induced DEGs ---------------------------------
SO <- readRDS(here("ALL_DCT.rds"))
Idents(SO) <- "analysis"

# DCT1
DCT1.kd <- FindMarkers(SO,
                       ident.1 = "DCT1_KD",
                       ident.2 = "DCT1_NK",
                       test.use = "DESeq2",
                       assay = "SCT",
                       min.pct = 0.1,
                       logfc.threshold = 0.3219)

# DCT2
DCT2.kd <- FindMarkers(SO,
                       ident.1 = "DCT2_KD",
                       ident.2 = "DCT2_NK",
                       test.use = "DESeq2",
                       assay = "SCT",
                       min.pct = 0.1,
                       logfc.threshold = 0.3219)

# Proliferation population
Prolif.kd <- FindMarkers(SO,
                         ident.1 = "Prolif_KD",
                         ident.2 = "Prolif_NK",
                         test.use = "DESeq2",
                         assay = "SCT",
                         min.pct = 0.1,
                         logfc.threshold = 0.3219)

# 4. DEGs of DCT1 and DCT2 subpopulation in control samples ----------------------------------

# Load CONTROL_DCT1 seurat object
SO <- readRDS(here("CONTROL_DCT1.rds"))
Idents(SO) <- "subclass.l2"

DCT1a_NK <- FindMarkers(SO,
                        ident.1 = "DCT1-A",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

DCT1b_NK <- FindMarkers(SO,
                        ident.1 = "DCT1-B",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

DCT1c_NK <- FindMarkers(SO,
                        ident.1 = "DCT1-C",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

DCT1d_NK <- FindMarkers(SO,
                        ident.1 = "DCT1-D",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

# Load CONTROL_DCT2 seurat object
SO <- readRDS(here("CONTROL_DCT2.rds"))
Idents(SO) <- "subclass.l2"

DCT2alpha.vs.beta_NK <- FindMarkers(SO,
                        ident.1 = "DCT2-alpha",
                        verbose = TRUE,
                        test.use = "DESeq2",
                        assay = "SCT",
                        min.pct = 0.1,
                        logfc.threshold = 0.3219)

# 5. Save all DEG lists -----------------------------------------------------
# data.frame
save(DCT1vsDCT2_NK, Prolif_NK, DCT1.kd, DCT2.kd, Prolif.kd, DCT1a_NK, DCT1b_NK, DCT1c_NK, DCT1d_NK, DCT2alpha.vs.beta_NK, file = here("ALL DEG files.RData"))
save(DCT1vsDCT2_NK, Prolif_NK, DCT1.kd, DCT2.kd, Prolif.kd, DCT1a_NK, DCT1b_NK, DCT1c_NK, DCT1d_NK, DCT2alpha.vs.beta_NK, file = here("ALL DEG files.rds"))

# excel file
dataset_names <- list('DCT1vsDCT2_NK' = DCT1vsDCT2_NK, 'Prolif_NK' = Prolif_NK, 'DCT1.kd' = DCT1.kd, 'DCT2.kd' = DCT2.kd, 'Prolif.kd' = Prolif.kd, 'DCT1a_NK' = DCT1a_NK, 'DCT1b_NK' = DCT1b_NK, 'DCT1c_NK' = DCT1c_NK, 'DCT1d_NK' = DCT1d_NK, 'DCT2alpha.vs.beta_NK' = DCT2alpha.vs.beta_NK)
write.xlsx(dataset_names, file = 'ALL DEGs.xlsx', rowNames=TRUE)
