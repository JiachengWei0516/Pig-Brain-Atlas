rm(list = ls())
library(ggdendro)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggtree)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("D:/Desktop")
dataDir ="D:/Desktop"
inputdir = "D:/Desktop"
outputdir = "D:/Desktop"

## Read the dataset
hipec <- readRDS(file = "mnn_figure1.rds")
hipec <- ScaleData(hipec, verbose = FALSE)

hipec <- FindVariableFeatures(hipec, selection.method = "vst", nfeatures = 3000)
hipec <- RunPCA(hipec, npcs = 30, verbose = FALSE)
group.by <- "celltype"


cls_order = c('Oligodendrocyte progenitor cells','Oligodendrocytes','Neurons','Astrocytes','Bergmann glia','Microglia','Perivascular macrophages','T cells','Vascular endothelial cells','Pericytes', 
              'Vascular smooth muscle cells','Ependymal cells','Pituitary cells','Olfactory ensheathing cells')

meta_data = hipec@meta.data
meta_data$group=meta_data$celltype
unique(meta_data$group)

## get the average expression of the dataset
if (TRUE){
  Idents(hipec) <- group.by
  avgs <- log(AverageExpression(hipec)$RNA + 1)
  saveRDS(avgs, file ="Avg.HIPEC.subtypes.rds") 
}
avgs <- readRDS(file ="Avg.HIPEC.subtypes.rds")

hvg <- rownames(hipec@reductions$pca@feature.loadings)  
avg_use <- avgs[hvg, ]


# Rectangular lines 
dendro <- dist(MinMax(scale(t(avg_use)), min = -4, max = 4)) %>% 
  hclust(., method = "ward.D2") %>% 
  as.dendrogram()


dendro <- reorder(dendro, match(colnames(avgs), rev(cls_order)), agglo.FUN=mean)
pdf("cluster_dendrogram.pdf", width = 7, height = 12)
ggdendrogram(dendro, rotate = TRUE)
dev.off()

theme_empty <- function(x) {
  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), legend.position = "none")
}

ddata <- dendro_data(dendro, type = "rectangle")
den_p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_empty()
den_p
ggsave("cluster.dendrogram_no_label.pdf",den_p,height = 7, width = 3)

theme_empty <- function(x) {
  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), legend.position = "none")
}

## Plot the region contribution
number_data  <- meta_data %>% 
  group_by(!!sym(group.by), group) %>% 
  summarize(ncells=n()) %>% 
  mutate(ncells2 = sqrt(ncells),
         ncells_div_1000 = ncells / 1000) %>% 
  ungroup()


Organ_data <- meta_data %>%
  group_by(!!sym(group.by), Organ) %>%
  summarize(ncells = n()) %>%
  ungroup() %>%
  group_by(Organ) %>%
  mutate(ratio = ncells/sum(ncells)) %>%
  ungroup() %>%
  group_by(!!sym(group.by)) %>%
  mutate(rratio = ratio * 100/sum(ratio))

Organ_data$Organ = factor(Organ_data$Organ, levels = c(
  'Olfactory bulb',
  'Cerebral cortex','Cingulate gyrus','Hippocampus',
  'Corpus callosum','Striatum',
  'Thalamus','Hypothalamus','Brain stem',
  'Pituitary',
  'Cerebellum',
  'Spinal cord'
))

library(paletteer)
reg_colors <- paste0(c(
  "#ffa15b","#ffc895","#ffe1c6",
  "#d7aed6", #红色
  "#958ec3","#b4b5d9", "#d5d5e9",  #紫色
  "#41a6d2", "#83c3de",  #蓝色
  "#88dbaa", "#b2e8b6", "#dcf5da"  #绿色
), "") %>% 
  setNames(c(
    'Pituitary','Cerebellum','Spinal cord',
    'Olfactory bulb',
    'Cerebral cortex','Cingulate gyrus','Hippocampus',
    'Corpus callosum','Striatum',
    'Thalamus','Hypothalamus','Brain stem'
  ))

celltype = c('VLMC','VEC','PC','VSMC','Neu','OPC','Oligo','Astro','BG','OEC','Epen','Pit','Micro','PVM','TC')

reversed_cls_order = rev(celltype)

number_data <- number_data %>%
  mutate(celltype = case_when(
    celltype == "Astrocytes" ~ "Astro",
    celltype == "Bergmann glia" ~ "BG",
    celltype == "Ependymal cells" ~ "EC",
    celltype == "Pericytes" ~ "PC",
    celltype == "Microglia" ~ "Micro",
    celltype == "T cells" ~ "TC",
    celltype == "Perivascular macrophages" ~ "PVM",
    celltype == "Neurons" ~ "Neu",
    celltype == "Oligodendrocyte progenitor cells" ~ "OPC",
    celltype == "Oligodendrocytes" ~ "Oligo",
    celltype == "Olfactory ensheathing cells" ~ "OEC",
    celltype == "Vascular and leptomeningeal cells" ~ "VLMC",
    celltype == "Vascular endothelial cells" ~ "VEC",
    celltype == "Vascular smooth muscle cells" ~ "VSMC",
    celltype == "Pituitary cells" ~ "Pit",
    TRUE ~ celltype
  ))


reg_p <- ggplot(Organ_data, aes_string(x = "rratio", y = group.by, fill = "Organ")) +
  geom_bar(positio = "stack", color = NA, stat = "identity") +
  scale_fill_manual(values = reg_colors) +
  scale_y_discrete(limits = rev(celltype))+
  scale_x_continuous(position = "top") +
  theme_empty() +
  theme(axis.text.x = element_text(size = 7, angle = 45,hjust = 0), axis.ticks.x = element_line(size = 0.2), axis.line.x = element_line(size = 0.2))

genes <- c("COL1A2","COL1A1", "CLDN5", "RGS5", "MYH11", "SNAP25", "PDGFRA","MBP", "AQP4", "NPY", "CFAP126", "GH1","PRL", "CALCR", "CD163", "PTPRC")

Idents(hipec) <- group.by
exp_p <- DotPlot(hipec, features = genes, cols = c("lightgrey", "red"), dot.min = 0.05, dot.scale = 3) +
  scale_x_discrete(limits = genes, position = "top") +
  scale_y_discrete(limits = rev(celltype)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 8,face = "bold",color = "black"),legend.text = element_text(size = 6,face = "bold",color = "black"), axis.title = element_blank(), axis.line = element_line(size = 0.2,color = "black"), axis.ticks = element_line(size = 0.2,color = "black"), axis.text.x = element_text(size = 7, face = "bold",angle = 45,hjust = 0,color = "black"), axis.text.y = element_blank())
pdf("MF1.Cluster.combo_gai1.pdf", width = 10, height = 6)
plot_grid(den_p, reg_p,  exp_p, nrow = 1, ncol = 3, rel_widths = c(0.5, 0.9, 2.5), align = "h") %>% print()
dev.off()

