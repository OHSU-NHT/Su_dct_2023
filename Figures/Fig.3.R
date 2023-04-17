
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(monocle3)
library(SeuratWrappers)
library(Signac)
library(patchwork)
library(Matrix)

# Fig.3D ---------------------------------------------------------------------------------------
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
  "F3D.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3E ---------------------------------------------------------------------------------------
SO <- PrepSCTFindMarkers(SO, verbose = TRUE)

markers <- FindAllMarkers(SO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

DotPlot(SO,
        features = unique(top10$gene),
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
  scale_y_discrete(limits = c("Prolif", "DCT2", "DCT1"))

ggsave(
  "F3E.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 10,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3F ---------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("DCT1_score1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) &
  scale_x_reverse()

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

# Fig.3G ---------------------------------------------------------------------------------------
FeaturePlot(SO,
            features = c("DCT2_score1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu"))) &
  scale_x_reverse()

ggsave(
  "F3G.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3H ---------------------------------------------------------------------------------------
# Covert Seurat object to a CellDataset object
cds <- as.cell_data_set(SO)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")

# find all possible partitions
all_partitions <- unique(cds@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]

# set all partitions to 1
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions %in% all_partitions] <- "1"

# set all cells to one partition, otherwise it won't show the trajectory across partitions. Set use_partition to FALSE as below.
cds <- learn_graph(cds, use_partition = F)

# Specify root cells programmatically
get_earliest_principal_node <- function(cds, time_bin="DCT"){ # set the time_bin as 'DCT' so that all cells are from the same bin, no bias.
  cell_ids <- which(colData(cds)[, "class"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

#Passing the programatically selected root node to order_cells() via the root_pr_node argument yields:
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds = cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1,
           label_groups_by_cluster=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 0.3) + scale_color_viridis_c() +
  scale_x_reverse()

ggsave(
  "F3H.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.3I
gene_list <- c("Erbb4", "Egf", "Trpm7", "Fgf13", "Col5a2", "Umod", "Ptgfr", "Stk32b", "Rtl4", "Abca13") # Top 10 in DCT1
pt.matrix <- exprs(cds)[match(gene_list,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- gene_list

tiff("F3I.tiff",
     units="in",
     width=7,
     height=3,
     res=700,
     compression = 'lzw'
     )

Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows       = "ward.D2",
  clustering_method_columns    = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

dev.off()

# Fig.3J
gene_list <- c("Slc8a1", "Arl15", "Calb1", "Slc2a9", "Phactr1", "Gls", "S100g", "Kl", "Klk1", "Egfem1")
pt.matrix <- exprs(cds)[match(gene_list,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- gene_list

tiff("F3J.tiff",
     units="in",
     width=7,
     height=3,
     res=700,
     compression = 'lzw'
)

Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows       = "ward.D2",
  clustering_method_columns    = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

dev.off()
