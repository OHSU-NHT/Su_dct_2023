
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)

# Fig.1D ---------------------------------------------------------------------------------------
SO <- readRDS(here("CONTROL_DCT.rds"))
Idents(SO) <- "subclass.l1"

mycols.DCT <- c("#66c2a5", # DCT1
                "#fc8d62", # DCT2
                "#A91554") # Prolif

DimPlot(SO, 
        reduction = "umap", 
        cols = mycols.DCT) +
  scale_x_reverse()

ggsave(
  "F1D.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)
