
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(tibble)

# Fig.3A -----------------------------------------------------------------------------------------------
SO <- readRDS(here("CONTROL_DCT.rds"))
Idents(SO) <- "subclass.l1"

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
  "F3A.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3B -----------------------------------------------------------------------------------------------
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
  "F3B.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3C -----------------------------------------------------------------------------------------------
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
  "F3C.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3D -----------------------------------------------------------------------------------------------
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
  "F3D_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
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
  "F3D_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3F -----------------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("Mg_score1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) +
  scale_x_reverse() +
  ggtitle("Magnesium score")

ggsave(
  "F3F.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3G -----------------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("Ca_score1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) +
  scale_x_reverse() +
  ggtitle("Calcium score")

ggsave(
  "F3F.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

