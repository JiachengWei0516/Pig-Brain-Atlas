library(Seurat)
library(monocle)
library(stringr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggsci)
setwd("D:/Desktop")
sce = readRDS("micro_E1E2_monocle.rds")
DefaultAssay(object = sce) <- "RNA"
Idents(sce)=sce@meta.data$celltype
set.seed(1)
sce = subset(sce, downsample = 3000)

# ----------- Define seurat2monocle function
Seurat2Monocle2=function(sce, do.init=T){
  expr_matrix=GetAssayData(sce, assay = 'RNA', slot = 'counts')
  cell_data=sce@meta.data
  gene_data=data.frame(
    gene_short_name=row.names(sce), #must have this column
    row.names = row.names(sce)
  )
  pd <- new("AnnotatedDataFrame", data = cell_data)
  fd <- new("AnnotatedDataFrame", data = gene_data)
  cds <- newCellDataSet(expr_matrix, 
                        phenoData = pd, 
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily=negbinomial.size())
  if(do.init){
    cds <- estimateSizeFactors(cds) #add column: Size_Factor
    cds <- estimateDispersions(cds)  
    cds=detectGenes(cds, min_expr = 0.1) #add column: num_cells_expressed
  }
  print( head(pData(cds)) )
  #
  return(cds)
}

prefix <- "micro_E1E2_celltype1_test_"
OUTPUT <- "D:/Desktop/"
setwd(OUTPUT)
cds = Seurat2Monocle2(sce)
class(cds)

sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
express_genes <- VariableFeatures(sce)
cds <- setOrderingFilter(cds, express_genes) 
write.table(express_genes,file=paste0(OUTPUT,"train.monocle.DEG.xls"),col.names=T,row.names=F,sep="\t",quote=F)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree') 

cds <- orderCells(cds)

# Visualization
cbPallete1 <- c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', 
                '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8',
                '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', 
                '#f7b6d2', '#dbdb8d')

pdf("train.monocle.pseudotime.pdf", width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "Pseudotime",size=1,show_backbone=TRUE)
dev.off()

p1= plot_cell_trajectory(cds, color_by = "State",cell_size = 1) 
ggsave(paste0(OUTPUT,"State_",prefix,".pdf"),p1)
ggsave(paste0(OUTPUT,"State_",prefix,".png"),p1)

p2=plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1)
ggsave(paste0(OUTPUT,"StateFacet_",prefix,".pdf"),p2,width = 16)
ggsave(paste0(OUTPUT,"StateFacet_",prefix,".png"),p2,width = 16)

p1= plot_cell_trajectory(cds, color_by = "celltype1") 
ggsave(paste0(OUTPUT,"celltype_",prefix,".pdf"),p1)
ggsave(paste0(OUTPUT,"celltype_",prefix,".png"),p1)

p2=plot_cell_trajectory(cds, color_by = "celltype1") +facet_wrap(~celltype, nrow = 1)
ggsave(paste0(OUTPUT,"celltypeFacet_",prefix,".pdf"),p2,width = 16)
ggsave(paste0(OUTPUT,"celltypeFacet_",prefix,".png"),p2,width = 16)

saveRDS(cds, paste0(OUTPUT,"cds", ".rds"))

cds <- orderCells(cds, root_state = 1)
saveRDS(cds,"micro_guanfang_pseudotime.rds")

p9 <- plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) + scale_color_gradient(low='#440154FF', high = '#FDE725FF')
ggsave(paste0(prefix,"_pseudotime.png"),p9,width=4,height=4)
ggsave(paste0(prefix,"_pseudotime.pdf"),p9,width=4,height=4)
setwd("D:/Desktop/")


wd <- "D:/Desktop"
prefix <- "micro_root_state"
cds <- readRDS(paste0(wd,prefix,"_pseudotime.rds"))
Time_diff <- differentialGeneTest(cds, cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,paste0(wd,prefix,"_Time_diff_all.csv"),row.names = F)
Time_diff <- read.csv('D:/Desktop/micro_state_gb2000_Time_diff_all.csv')
rownames(Time_diff) <- Time_diff$gene_short_name

sig_gene_names <- subset(Time_diff, qval < 0.01)
Time_genes <- sig_gene_names %>% pull(gene_short_name) %>% as.character()
p11=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=5, show_rownames=T, return_heatmap=T)


ggsave(paste0(wd,prefix,"_Time_heatmapAll_clut",clut,".png"),p11)
ggsave(paste0(wd,prefix,"_Time_heatmapAll_clut",clut,".pdf"),p11,width=5,height=10)

p11_clusters <- as.data.frame(cutree(p11$tree_row, k = 6))
names(p11_clusters) <- "cluster"
p11_clusters$gene <- rownames(p11_clusters)
tmp <- left_join(sig_gene_names,p11_clusters, by=c("gene_short_name"="gene"))
write.csv(tmp, paste0(wd,prefix,"_Time_heatmapAll_clut_sig.csv"), row.names = F)

saveRDS(cds,paste0(wd,prefix,"_pseudotime_clut.rds"))
