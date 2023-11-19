
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)

# Fig.6A-D. DCT1 subpopulation ------------------------------------------------------------------------------------------------
# Fig.6A ----------------------------------------------------------------------------------------------------------------------
SO <- readRDS(here("CONTROL_DCT1.rds"))
Idents(SO) <- "subclass.l2"
my_levels <- c("DCT1-A", "DCT1-B", "DCT1-C", "DCT1-D")
Idents(SO) <- factor(x = Idents(SO), levels = my_levels)

mycols.DCT1.sub <- c("#56ba5a", # DCT1-A
                     "#52b7bd", # DCT1-B
                     "#A8D8B5", # DCT1-C
                     "#4d8185") # DCT1-D
                     
DimPlot(SO,
        reduction = "umap",
        label = TRUE,
        pt.size = 0,
        label.size = 6,
        cols = mycols.DCT1.sub) +
  NoLegend()

ggsave(
  "F5A.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6B ----------------------------------------------------------------------------------------------------------------------
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
  scale_y_discrete(limits = c("DCT1-D", "DCT1-C", "DCT1-B", "DCT1-A"))

ggsave(
  "F5B.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 6.8,
  height = 3.3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6C ----------------------------------------------------------------------------------------------------------------------
Mg.markers <- c("Slc41a3", "Prox1", "Umod", "Cnnm2", "Fxyd2", "Egf", "Trpm6", "Trpm7")
DotPlot(SO,
        features = Mg.markers,
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
  coord_flip() + 
  RotatedAxis() +
  ggtitle("Magnesium") 

ggsave(
  "F5C_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

Ca.markers <- c ("Ryr2", "Trpv5", "S100g", "Vdr", "Calb1", "Slc8a1")
DotPlot(SO,
        features = Ca.markers,
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
  coord_flip() + 
  RotatedAxis() +
  ggtitle("Calcium")

ggsave(
  "F5C_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

Na.markers <- c("Hsd11b2", "Scnn1g", "Scnn1b", "Scnn1a", "Klk1", "Nr3c2", "Stk39", "Klhl3", "Wnk4", "Wnk1", "Slc12a3")
DotPlot(SO,
        features = Na.markers,
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
  coord_flip() + 
  RotatedAxis() +
  ggtitle("Sodium")

ggsave(
  "F5C_3.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6D ----------------------------------------------------------------------------------------------------------------------
DotPlot(SO,
        features = c("Hoxb7", "Hoxd10", "Emx1", "Sall3"),
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  coord_flip() + 
  RotatedAxis()

ggsave(
  "F5D_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 2,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

DotPlot(SO,
        features = c("Hoxb7", "Hoxd10"),
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  coord_flip() + 
  RotatedAxis()

ggsave(
  "F5D_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 1.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6E-H. DCT2 subpopulation ------------------------------------------------------------------------------------------------
# Fig.6E ----------------------------------------------------------------------------------------------------------------------
SO <- readRDS(here("CONTROL_DCT2.rds"))
Idents(SO) <- "subclass.l2"


mycols.DCT2.sub <- c("#CC79A7", # DCT2-X
                     "#D55E00") # DCT2-Y
                              
DimPlot(SO,
        reduction = "umap",
        label = TRUE,
        pt.size = 0,
        label.size = 6,
        cols = mycols.DCT2.sub) +
  NoLegend()

ggsave(
  "F5E.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6F ----------------------------------------------------------------------------------------------------------------------
SO <- PrepSCTFindMarkers(SO, verbose = TRUE)
ALL.markers.DCT2.sub <- FindAllMarkers(SO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ALL.markers.DCT2.sub %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- ALL.markers.DCT2.sub %>% group_by(cluster) %>% top_n(5, avg_log2FC)

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
  scale_y_discrete(limits = c("DCT2-Y", "DCT2-X"))

ggsave(
  "F5F.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 6,
  height = 3.3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6G ----------------------------------------------------------------------------------------------------------------------
Mg.markers <- c("Slc41a3", "Prox1", "Umod", "Cnnm2", "Fxyd2", "Egf", "Trpm6", "Trpm7")
DotPlot(SO,
        features = Mg.markers,
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
  coord_flip() + 
  RotatedAxis() +
  ggtitle("Magnesium") 

ggsave(
  "F5G_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

Ca.markers <- c ("Ryr2", "Trpv5", "S100g", "Vdr", "Calb1", "Slc8a1")
DotPlot(SO,
        features = Ca.markers,
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
  coord_flip() + 
  RotatedAxis() +
  ggtitle("Calcium")

ggsave(
  "F5G_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

Na.markers <- c("Hsd11b2", "Scnn1g", "Scnn1b", "Scnn1a", "Klk1", "Nr3c2", "Stk39", "Klhl3", "Wnk4", "Wnk1", "Slc12a3")
DotPlot(SO,
        features = Na.markers,
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
  coord_flip() + 
  RotatedAxis() +
  ggtitle("Sodium")

ggsave(
  "F5G_3.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6H ----------------------------------------------------------------------------------------------------------------------
DotPlot(SO,
        features = c("Hoxb7", "Hoxd10", "Emx1", "Sall3"), 
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  coord_flip() + 
  RotatedAxis() + 
  theme(axis.title = element_blank())

ggsave(
  "F5H_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 2,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

DotPlot(SO,
        features = c("Hoxb7", "Hoxd10"), 
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  coord_flip() + 
  RotatedAxis() + 
  theme(axis.title = element_blank())

ggsave(
  "F5H_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 1.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)
