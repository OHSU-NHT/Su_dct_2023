
library(Seurat)
library(here)
library(monocle3)
library(SeuratWrappers)
library(Signac)
library(patchwork)
library(Matrix)
library(ggplot2)

# Building trajectories with Monocle 3
SO <- readRDS(here("CONTROL_DCT.rds"))

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
           cell_size = 0.1) + scale_color_viridis_c() +
  scale_x_reverse()

# Add Pseudotime MetaData to the Seurat object
SO <- AddMetaData(object = SO,
                  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
                  col.name = "Pseudotime")


