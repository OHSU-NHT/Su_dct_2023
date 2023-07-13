
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(tibble)
library(clusterProfiler)

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

# Fig.3E -----------------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("Mg_score1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) +
  scale_x_reverse() +
  ggtitle("Magnesium score")

ggsave(
  "F3E.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3F -----------------------------------------------------------------------------------------------
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

# Fig.3G-H Pathway analysis - Import the 'ALL DEG files.RData' file (supplemental data from original publication) into the R Studio environment
# Convert Gene Symbols to ENTREZ IDs
markers <- DCT1vsDCT2_NK %>% rownames_to_column(var="SYMBOL")

ENTREZ_list <- bitr(geneID = rownames(DCT1vsDCT2_NK),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db"
)

markers <-  ENTREZ_list %>% inner_join(markers, by = "SYMBOL")

# Remove genes that are not statistically significant 
markers <-  markers %>% dplyr::filter(p_val_adj < 0.05)
head(markers, n = 50)

# Fig.3G DCT1 pathways ----------------------------------------------------------------------------------

# Up-regulated Genes (DCT1-enriched) = markers with positive fold change values 
pos.markers <-  markers %>% dplyr::filter(avg_log2FC > 0) %>%  arrange(desc(abs(avg_log2FC)))
head(pos.markers, n = 50)
pos.ranks <- pos.markers$ENTREZID[abs(markers$avg_log2FC) > 0.3219]
head(pos.ranks)

# GO - Biological Processes 
pos_go.BP <- enrichGO(pos.ranks,
                      OrgDb = "org.Mm.eg.db",
                      ont ="BP",
                      readable=TRUE)

# Figure 3G - barplot of DCT1-enriched pathway
df.pos_go.BP <- data.frame(GO_Biological_Process = pos_go.BP@result$Description,
                           Count = pos_go.BP@result$Count,
                           Pvalue = pos_go.BP@result$p.adjust)
df.pos_go.BP <- df.pos_go.BP[order(df.pos_go.BP[,"Count"], decreasing = TRUE),]
df.pos_go.BP$observation <- 1:nrow(df.pos_go.BP)
df.pos_go.BP <- df.pos_go.BP %>% mutate(highlight = dplyr::case_when(observation <=36 ~ "Yes",
                                                                     TRUE ~ "No"))
top80 <- df.pos_go.BP %>% top_n(80, Count)

ggplot(top80, aes(x = reorder(GO_Biological_Process, observation, decreasing = FALSE),
                  y = Count,
                  fill = highlight
)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_manual(values=c("#999999", "royal blue")) +
  ylab("Gene Count") +
  xlab("GO Biological Process") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 1.5,
                                   size = 15),
        axis.text.y = element_text(size = 15))

ggsave(
  "F3G.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 10,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3H DCT2 pathways ----------------------------------------------------------------------------------

# Down-regulated Genes (DCT2-enriched) = markers with positive fold change values 
neg.markers <-  markers %>% dplyr::filter(avg_log2FC < 0) %>%  arrange(desc(abs(avg_log2FC)))
head(neg.markers, n = 50)
neg.ranks <- neg.markers$ENTREZID[abs(markers$avg_log2FC) > 0.3219]
head(neg.ranks)

# GO - Biological Processes 
neg_go.BP <- enrichGO(neg.ranks,
                      OrgDb = "org.Mm.eg.db",
                      ont ="BP",
                      readable=TRUE)

# Figure 3H - barplot of DCT2-enriched pathway
df.neg_go.BP <- data.frame(GO_Biological_Process = neg_go.BP@result$Description,
                           Count = neg_go.BP@result$Count,
                           Pvalue = neg_go.BP@result$p.adjust)
df.neg_go.BP <- df.neg_go.BP[order(df.neg_go.BP[,"Count"], decreasing = TRUE),]
df.neg_go.BP$observation <- 1:nrow(df.neg_go.BP)
df.neg_go.BP <- df.neg_go.BP %>% mutate(highlight = dplyr::case_when(observation <= 30 ~ "Yes",
                                                                     TRUE ~ "No"))
top80 <- df.neg_go.BP %>% top_n(80, Count)

ggplot(top80, aes(x = reorder(GO_Biological_Process, observation, decreasing = FALSE),
                  y = Count,
                  fill = highlight
)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_manual(values=c("#999999", "royal blue")) +
  ylab("Gene Count") +
  xlab("GO Biological Process") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 1.5,
                                   size = 15),
        axis.text.y = element_text(size = 15))

ggsave(
  "F3H.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 11,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)
