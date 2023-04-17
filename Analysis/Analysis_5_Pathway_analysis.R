
library(Seurat)
library(dplyr)
library(tibble)
library(clusterProfiler)

# Import the 'ALL DEG files.RData' file (supplemental data from original publication) into the R Studio environment

# 1. Convert Gene Symbols to ENTREZ IDs ------------------------------------------------------------------
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

# 2. DCT1 pathways ---------------------------------------------------------------------------

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

# Barplot of DCT1-enriched pathway
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

# 3. DCT2 pathways ---------------------------------------------------------------------------

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

# Barplot of DCT2-enriched pathway
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
