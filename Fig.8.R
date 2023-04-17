
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)

# Fig. 8A ------------------------------------------------------------------------------------------------
SO <- readRDS(here("ALL_DCT.rds"))
Idents(SO) <- "subclass.l1"
my_levels <- c("DCT1", "DCT2", "Prolif")
Idents(SO) <- factor(x = Idents(SO), levels = my_levels)

mycols.DCT <- c("#66c2a5", # DCT1
                "#fc8d62", # DCT2
                "#A91554") # Prolif
                         
DimPlot(SO,
        reduction = "umap",
        pt.size = 0,
        cols = mycols.DCT,
        split.by = "Diet") &
  scale_x_reverse()

ggsave(
  "F8A.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 6,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 8G ------------------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("DCT_score1")) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) &
  scale_x_reverse()

ggsave(
  "F8G_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

FeaturePlot(SO,
            features = c("DCT_score1"),
            split.by = "Diet") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) &
  scale_x_reverse() &
  NoAxes()

ggsave(
  "F8G_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 8,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)


