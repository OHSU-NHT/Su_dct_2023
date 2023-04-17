
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(tibble)
library(ggpubr)

# Fig. 9A ------------------------------------------------------------------------------------------------
SO <- readRDS(here("ALL_DCT.rds"))
Idents(SO) <- "subclass.l1"
my_levels <- c("DCT1", "DCT2", "Prolif")
Idents(SO) <- factor(x = Idents(SO), levels = my_levels)

VlnPlot(SO,
        features = c("Ca_score1"),
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
  "F9A.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig. 9B ------------------------------------------------------------------------------------------------
df <- FetchData(object = SO, vars = c("subclass.l1", "analysis", "Rep", "Diet", "Ca_score1")) %>% rownames_to_column()
df2 <- subset(df, subclass.l1 == c("DCT1", "DCT2"))
df2$analysis <- factor(x = df2$analysis, levels = c("DCT1_NK", "DCT2_NK", "DCT1_KD", "DCT2_KD"))
df$Diet <- factor(x = df$Diet, levels = c("NK", "KD"))

ggplot(df2, aes(x = Ca_score1, group = analysis, fill = analysis)) + 
  geom_density(adjust=1.5,
               alpha=.6,
               color = c("black")) +
  scale_fill_manual(values = c("white", "lightgrey", "lightskyblue", "#4CC3FF"))+
  theme_pubr(legend = "right") +
  ylab("Probability Density") +
  xlab("Calcium Score") +
  coord_flip()

ggsave(
  "F9B.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)
