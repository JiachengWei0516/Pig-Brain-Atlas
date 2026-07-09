
library(Seurat)
library(tidyverse)
library(patchwork)
setwd('D:/Desktop')

rm(list=ls())

## pig data
pig1=readRDS('pig.rds')
Idents(pig1)=pig1@meta.data$celltype_pig
saveRDS(pig1,'pig1.rds')


## human
human1 = readRDS("human.rds")
Idents(human1)=human1@meta.data$celltype_human

## mouse
mouse1 = readRDS('mouse.rds')
Idents(mouse1)=mouse1@meta.data$celltype_mouse

commonGenes <- intersect(rownames(pig1),rownames(human1))
commonGenes = intersect(commonGenes, rownames(mouse1))

dat_pig <- subset(pig1, features = as.matrix(commonGenes))
dat_Human_PFC <- subset(human1, features = as.matrix(commonGenes))
dat_mouse <- subset(mouse1, features = as.matrix(commonGenes))

dat_pig <- FindVariableFeatures(dat_pig, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
dat_Human_PFC <- FindVariableFeatures(dat_Human_PFC, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
dat_mouse <- FindVariableFeatures(dat_mouse, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

endolist <- NULL
endolist[[1]] <- dat_pig
endolist[[2]] <- dat_human
endolist[[3]] = dat_mouse

##1.3 Integreated Data 
options(future.globals.maxSize=1048576000000)
features <- SelectIntegrationFeatures(object.list=endolist)
endoAnchors <- FindIntegrationAnchors(object.list=endolist,anchor.features=features)
endoCombined <- IntegrateData(anchorset = endoAnchors)
Instcell.combined <- endoCombined
DefaultAssay(Instcell.combined) <- "integrated"

Instcell.combined <- ScaleData(Instcell.combined, verbose = FALSE)
Instcell.combined <- RunPCA(Instcell.combined, npcs = 30, verbose = FALSE)
Instcell.combined <- RunUMAP(Instcell.combined, reduction = "pca", dims = 1:30)
Instcell.combined <- FindNeighbors(Instcell.combined, reduction = "pca", dims = 1:30)
Instcell.combined <- FindClusters(Instcell.combined, resolution = 1.5)
saveRDS(Instcell.combined,file='pig_Human_Hippo_Instcell_combined_res_1.5_final.rds')

mnn_batch <- data.frame(mnn_batch =c(rep("pig",8851),rep("human",9017),rep("mouse",9564))) 

human = subset(x=Instcell.combined,subset= (mnn_batch == "human") )

mycolour = c("#84b3d4","#fccde5","#d43230","#b25c2c","#33a02c","#6d419c","#d1bcdb")
p4 <- DimPlot(human, reduction = "umap",cols=mycolour, group.by = "celltype_merge")
ggsave("inter_human_white.pdf", plot = p4, width = 5, height = 4)

mouse = subset(x=Instcell.combined,subset= (mnn_batch == "mouse") )
p5 <- DimPlot(mouse, reduction = "umap",cols=mycolour, group.by = "celltype_merge",na.value = "#00000000")
ggsave("inter_mouse_white.pdf", plot = p5,width = 5, height = 4)

pig = subset(x=Instcell.combined,subset= (mnn_batch == "pig") )
p6 <- DimPlot(pig, reduction = "umap",cols=mycolour, group.by = "celltype_merge",na.value = "#00000000")
ggsave("inter_pig_white.pdf", plot = p6,width = 5, height = 4)


mycolour = c("pig" = "#a6cee3",
             "Human_PFC" = "#ffd92f",
             "mouse" = "#fb9a99")
p7 <- DimPlot(Instcell.combined, reduction = "umap",cols=mycolour, group.by = "mnn_batch")
ggsave("pig_human_mouse_mnn_batch.pdf", plot = p7,width = 5, height = 4)


rm(list=ls())
library(Seurat)
library(dplyr)
library(circlize)
library(reshape2)
library(ComplexHeatmap)
setwd('D:/Desktop/')

Instcell.combined <- readRDS('pig_Human_Hippo_Instcell_combined_res_1.5_final.rds')

Instcell.combined$species <- Instcell.combined$mnn_batch
Instcell.combined$species_celltype <- paste0(Instcell.combined$species, "_", Instcell.combined$celltype_merge)

# =====================
# 设置颜色映射
# =====================
celltype_color_map <- c(
  "Astro" = "#84b3d4",
  "Micro"  = "#33a02c",
  "Oligo"      = "#6d419c",
  "OPC"        = "#d1bcdb",
  "DG GC"    = "#fccde5",
  "ExN"      = "#d43230",
  "InN"        = "#b25c2c"
)

link_color_map <- function(corr) {
  ifelse(corr > 0.9, "#7f0000",           
         ifelse(corr > 0.8, "#fdd49e",    
                ifelse(corr > 0.5, "#dde8f2", "#d8d8d8"))) 
}

avg_exp_list <- AverageExpression(Instcell.combined, group.by = "species_celltype", return.seurat = FALSE)
avg_exp_mat <- as.matrix(avg_exp_list$integrated)  # 确保是矩阵格式

cor_matrix <- cor(avg_exp_mat, method = "pearson")
write.csv(cor_matrix,"RNA_cor_matrix_finall.csv")

cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("from", "to", "correlation")


df_long <- cor_df %>%
  filter(from != to) %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(from, to)), collapse = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(from, to, correlation)

extract_species <- function(x) strsplit(x, "-")[[1]][1]


df_long$species_from <- sapply(as.character(df_long$from), extract_species)
df_long$species_to   <- sapply(as.character(df_long$to), extract_species)

df_long <- df_long %>% filter(species_from != species_to)

df_long$link_color <- link_color_map(df_long$correlation)

all_sectors <- unique(c(df_long$from, df_long$to))
sector_colors <- setNames(sapply(all_sectors, function(x) {
  x <- as.character(x)  #
  if (is.na(x)) {
    warning("⚠️")
    return("#999999")
  }
  parts <- unlist(strsplit(x, "-"))
  ct <- tail(parts, 1)  
  if (ct %in% names(celltype_color_map)) {
    celltype_color_map[[ct]]
  } else {
    warning(paste("❗ 未知的细胞类型:", ct, "请检查 all_sectors 中的命名"))
    "#999999"
  }
}), all_sectors)

pdf("cross_species_chord_diagram_new_1.pdf", width = 4, height = 4)
circos.clear()
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(
  x = df_long[, c("from", "to", "correlation")],
  grid.col = sector_colors,
  col = df_long$link_color,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.05),
  transparency = 0.4,
  directional = 0
)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1,col='black')
}, bg.border = NA)

lgd1 <- Legend(title = "Cell Type", labels = names(celltype_color_map),
               legend_gp = gpar(fill = celltype_color_map))

lgd2 <- Legend(title = "Correlation", 
               labels = c(">0.9", "0.8–0.9", "0.5–0.8", "<0.5"),
               legend_gp = gpar(fill = c("#7f0000", "#fdd49e", "#dde8f2", "#d8d8d8")))

draw(packLegend(lgd1, lgd2), x = unit(0.8, "npc"), y = unit(0.9, "npc"), just = c("left", "top"))


all_sectors <- get.all.sector.index()

extract_species <- function(x) {
  strsplit(x, "-")[[1]][1]
}

sector_species_map <- setNames(sapply(all_sectors, extract_species), all_sectors)
species_sectors <- split(names(sector_species_map), sector_species_map)

highlight_colors <- c(
  "pig" = "#a6cee3",
  "Human" = "#ffd92f",
  "mouse" = "#fb9a99"
)

for (sp in names(species_sectors)) {
  highlight.sector(
    sector.index = species_sectors[[sp]],
    track.index = 1,
    col = highlight_colors[sp],
    text = sp,
    cex = 1,
    niceFacing = TRUE,
    text.col = "black",
    text.vjust = 0.5
  )
}
dev.off()
