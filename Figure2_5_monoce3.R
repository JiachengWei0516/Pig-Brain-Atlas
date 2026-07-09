library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(reticulate)
library(SeuratWrappers)

rm(list=ls())

setwd('D:/Desktop/')

prefix <- "mnn_oligo"	
scRNAsub <- readRDS('mnn_oligo_main_harmony.rds')

#导入umap信息
sc_embedding <- read.table("sc_embedding.csv", header = TRUE, sep = "\t", row.names = 1)
scRNAsub[["umap_sc"]] <- CreateDimReducObject(as.matrix(sc_embedding),key = "umap_")
saveRDS(scRNAsub,'oligo_monocle_umap_meta.rds')

data <- GetAssayData(scRNAsub, assay = 'RNA', layer = 'counts')
cell_metadata <- scRNAsub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 30)

cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNAsub, reduction = "Xumap_")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')

ggsave(paste0(prefix,"_cds_umap.png"), plot = p1, width = 8, height = 7)
ggsave(paste0(prefix,"_seurat_umap.png"), plot = p2, width = 8, height = 7)

cds <- cluster_cells(cds,cluster_method = "leiden",resolution =0.00001) 
p3 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p4 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p3, p4)
ggsave(paste0(prefix,"_partitionID.png"), plot = p, width = 12, height = 7)


cds <- learn_graph(cds)
p5 = plot_cells(cds, color_cells_by="celltype",label_groups_by_cluster = FALSE, label_leaves = TRUE, 
                label_branch_points = TRUE,graph_label_size=3)
ggsave(paste0(prefix,"_graph_0.00005.png"), plot = p5, width = 8, height = 7)
saveRDS(cds,paste0(prefix,"_learn_graph.rds"))

sub = subset(scRNAsub, celltype == 'OPC')
embed <- data.frame(Embeddings(sub, reduction = "Xumap_"))
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
p6 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  
                 label_branch_points = FALSE,label_groups_by_cluster =FALSE,label_roots =FALSE,show_trajectory_graph =FALSE,
                 cell_size = 0.1)
ggsave(paste0(prefix,"_pseudotime_root.pdf"),plot=p6,width=5,height=4)


cds <- order_cells(cds)
p6 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE,label_groups_by_cluster =FALSE,label_roots =FALSE,show_trajectory_graph =FALSE,
                 cell_size = 0.1)
p6$layers[[1]]$aes_params$stroke <- 0

ggsave("oligo_pseudotime_root.pdf",plot=p6,width=5,height=4)