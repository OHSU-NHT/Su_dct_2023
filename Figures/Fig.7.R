
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(tibble)

# Fig. 7C ------------------------------------------------------------------------------------------------
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
        cols = mycols.DCT) +
  scale_x_reverse()

ggsave(
  "F7C_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

DimPlot(SO,
        reduction = "umap",
        pt.size = 0,
        group.by = "Rep") +
  scale_x_reverse()

ggsave(
  "F7C_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 7D ------------------------------------------------------------------------------------------------
VlnPlot(SO,
        features = c("DCT1_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F7D_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

VlnPlot(SO,
        features = c("DCT2_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F7D_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

VlnPlot(SO,
        features = c("ENaC_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F7D_3.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 7E ------------------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("DCT1_score1"),
            split.by = "Diet") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) &
  scale_x_reverse() &
  NoAxes()

ggsave(
  "F7E_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 7.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

FeaturePlot(SO,
            features = c("DCT2_score1"),
            split.by = "Diet") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) &
  scale_x_reverse() &
  NoAxes()

ggsave(
  "F7E_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 7.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 7F-H ----------------------------------------------------------------------------------------------
df <- FetchData(object = SO, vars = c("subclass.l1", "analysis", "Rep", "Diet", "DCT1_score1", "DCT2_score1", "ENaC_score1")) %>% rownames_to_column()
df2 <- subset(df, subclass.l1 == c("DCT1", "DCT2"))
df2$analysis <- factor(x = df2$analysis, levels = c("DCT1_NK", "DCT2_NK", "DCT1_KD", "DCT2_KD"))
df$Diet <- factor(x = df$Diet, levels = c("NK", "KD"))

# Fig. 7F ------------------------------------------------------------------------------------------------
ggplot(df2, aes(x = DCT1_score1, group = analysis, fill = analysis)) + 
  geom_density(adjust=1.5,
               alpha=.6,
               color = c("black")) +
  scale_fill_manual(values = c("white", "lightgrey", "lightskyblue", "#4CC3FF"))+
  theme_pubr(legend = "right") +
  ylab("Probability Density") +
  xlab("DCT1 Score") +
  coord_flip()

ggsave(
  "F7F.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 7G ------------------------------------------------------------------------------------------------
ggplot(df2, aes(x = DCT2_score1, group = analysis, fill = analysis)) + 
  geom_density(adjust=1.5,
               alpha=.6,
               color = c("black")) +
  scale_fill_manual(values = c("white", "lightgrey", "lightskyblue", "#4CC3FF"))+
  theme_pubr(legend = "right") +
  ylab("Probability Density") +
  xlab("DCT2 Score") +
  coord_flip()

ggsave(
  "F7G.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 7H ------------------------------------------------------------------------------------------------
ggplot(df2, aes(x = ENaC_score1, group = analysis, fill = analysis)) + 
  geom_density(adjust=1.5,
               alpha=.6,
               color = c("black")) +
  scale_fill_manual(values = c("white", "lightgrey", "lightskyblue", "#4CC3FF"))+
  theme_pubr(legend = "right") +
  ylab("Probability Density") +
  xlab("ENaC Score") +
  coord_flip()

ggsave(
  "F7H.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)
