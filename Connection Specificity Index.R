library(Seurat)
library(ggplot2)
library(dplyr)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(pheatmap)
library(circlize)
library(gridExtra)
library(reshape2)


pbmc3k_sub = readRDS('/home/snatac/02.wjc/09.scanpy/mnn_counts_astro.rds')

loom <- open_loom("/home/snatac/02.wjc/09.scanpy/00.pyscenic/out_SCENIC_Astro_mc9nr.loom")

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)

regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')

regulonAucThresholds <- get_regulon_thresholds(loom)
sub_regulonAUC <- regulonAUC[,match(colnames(pbmc3k_sub),colnames(regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(pbmc3k_sub))

pbmc3k_sub@meta.data = cbind(pbmc3k_sub@meta.data ,t(sub_regulonAUC@assays@data$AUC))

rss<-calcRSS(AUC=getAUC(sub_regulonAUC),cellAnnotation=pbmc3k_sub@meta.data$celltype8)
pdf('rss_plot_Micro_state.pdf')
rssPlot <- plotRSS(rss)
rssPlot$plot
dev.off()

auc <- getAUC(sub_regulonAUC)

auc_mat <- as.matrix(auc)  # regulon x cell
if (exists("reg_thresh")) {

  reg_thresh <- reg_thresh[rownames(auc_mat)]
  auc_bin <- (auc_mat > reg_thresh) * 1
} else {
  auc_z <- t(scale(t(auc_mat)))        
  auc_bin <- (auc_z > 0) * 1
}

# TF names:
tf_names <- rownames(auc_bin)
nTF <- length(tf_names)

pcc_mat <- cor(t(auc_bin), method = "pearson")

csi_mat <- matrix(0, nTF, nTF, dimnames = list(tf_names, tf_names))

for (i in seq_len(nTF)) {
  if (i %% 50 == 0) message("progress:", i, "/", nTF)
  for (j in seq_len(nTF)) {
    if (i == j) { csi_mat[i,j] <- 1; next }
    p_ij <- pcc_mat[i,j]
    others <- setdiff(seq_len(nTF), c(i,j))
    cond <- (pcc_mat[others, i] < p_ij) & (pcc_mat[others, j] < p_ij)
    csi_mat[i,j] <- sum(cond, na.rm = TRUE) / length(others)
  }
}

csi_sym <- (csi_mat + t(csi_mat)) / 2
csi_mat <- csi_sym

p <- pheatmap(csi_mat,
              cluster_rows = TRUE, cluster_cols = TRUE,
              show_rownames = FALSE, show_colnames = FALSE,

              main = "TF-TF CSI heatmap")
pdf('TF-TF_CSI_heatmap.pdf',width=8,height=8)  
p
dev.off()
clustered_tf_index <- p$tree_row$order

heatmap_tf_order <- rownames(csi_mat)[clustered_tf_index]

print(heatmap_tf_order)

write.table(data.frame(TF_name = heatmap_tf_order, 
                       Heatmap_Row_Number = 1:length(heatmap_tf_order)),
            file = "heatmap_TF_order.txt", 
            row.names = FALSE, sep = "\t")



# ---------- module activity----------
modules <- list(
  M1 = c("AR(+)", "IRF4(+)", "TP73(+)", "FOXA1(+)", "MYB(+)", "NKX3-1(+)", 
         "SP4(+)", "EGR2(+)", "ZBTB48(+)"),
  M2 = c("E2F1(+)", "RFX3(+)", "SOX6(+)", "BCL11A(+)", "SP8(+)"),
  M3 = c("TET1(+)", "ZEB1(+)", "RAD21(+)", "SP3(+)", "NR2C2(+)",
         "POU2F1(+)", "DMRTA2(+)", "NR1H3(+)"),
  M4 = c("KLF5(+)", "PAX3(+)", "MEIS1(+)", "EN2(+)", "SOX18(+)", "ELF5(+)", "ETV1(+)"),
  M5 = c("JUNB(+)", "NEUROD2(+)", "PPARD(+)", "TRAF4(+)", "ONECUT1(+)", 
         "LHX1(+)", "UNCX(+)", "BCLAF1(+)", "TAF7(+)", "TBX6(+)", 
         "MAZ(+)", "PPARG(+)", "CREB3L4(+)", "POU3F4(+)"),
  M6 = c("LHX5(+)", "MSX2(+)", "CREB1(+)", "ZNF524(+)", "BCL3(+)", "ATF1(+)", "E2F4(+)"),
  M7 = c("SIX3(+)", "DLX1(+)", "EN1(+)", "ZNF471(+)", "BRF1(+)", "IRF3(+)", "HOXA2(+)", "IRX5(+)"),  
  M8 = c("NR4A2(+)", "SOX13(+)", "ATF4(+)", "EZH2(+)", "E2F8(+)", "IRF1(+)"),
  M9 = c("ZNF768(+)", "MAFG(+)", "PPARA(+)", "E2F3(+)", "GABPB1(+)", "TBX3(+)", "BCL6B(+)", "CEBPD(+)", "CHD1(+)",
         "ZNF263(+)", "ZNF274(+)", "NFYA(+)", "ZNF143(+)", "BDP1(+)", "ERF(+)", "DDIT3(+)", "DMRTA1(+)", "FOXK1(+)",
         "GATA5(+)", "ZNF467(+)", "GMEB2(+)", "NFIL3(+)"),
  M10 = c("ELF1(+)", "MLXIPL(+)", "ELK3(+)", "ELF4(+)", "MAF(+)", "ETS2(+)", "IRF8(+)",
          "PRDM1(+)", "CEBPA(+)", "CEBPB(+)", "IRF5(+)", "IKZF1(+)", "MEF2A(+)",
          "ZNF683(+)", "NFKB1(+)", "RELA(+)", "EHF(+)", "HESX1(+)", "POU2F2(+)",
          "ETV6(+)", "SOX17(+)", "FOSL2(+)", "RUNX1(+)", "MAFF(+)", "ATF3(+)",
          "FOS(+)","LHX9(+)", "IRF6(+)", "MEF2D(+)"),
  M11 = c("RXRG(+)", "NFE2L1(+)", "TCF7(+)",
          "MITF(+)", "OTX1(+)", "FOXO1(+)", "STAT1(+)", "HLF(+)", "FOXN2(+)", "LHX2(+)", "SOX9(+)", "ETV4(+)", "EMX2(+)",
          "SREBF1(+)", "TCF12(+)", "NR3C2(+)", "RORA(+)", "YBX1(+)", "NF1(+)", "ZBTB33(+)", "KAT2A(+)", "STAT2(+)", "NR2F6(+)",
          "SREBF2(+)", "VEZF1(+)", "ZMIZ1(+)", "TBL1XR1(+)", "RBBP9(+)", "ZBTB41(+)"),
  M12 = c("OTX2(+)", "SHOX2(+)", "HEYL(+)", "ZNF335(+)", "BCL6(+)", 
          "ZBTB2(+)")
)

module_activity <- sapply(modules, function(tf_list) {
  
  tf_list <- intersect(tf_list, rownames(auc))
  if (length(tf_list) == 0) return(rep(NA, ncol(auc)))
  
  colMeans(auc[tf_list, , drop = FALSE])
})

module_activity <- as.data.frame(module_activity)
pbmc3k_sub = readRDS('/home/snatac/02.wjc/09.scanpy/mnn_counts_astro.rds')
pbmc3k_sub@meta.data <- cbind(pbmc3k_sub@meta.data, module_activity)

pbmc3k_sub <- readRDS('D:/Desktop/Astro_module_activity_M1_M12.rds')
my_colors <- c("#f0efef", "#cecece", "#d6a933", "#3a1b20")

my_breaks <- seq(0, 0.25, length.out = 6)
FeaturePlot(
  pbmc3k_sub,
  features = "M2",
  reduction = "Xumap_"
) +
  scale_color_gradientn(
    colours = my_colors,
    limits = c(0, 0.25),         # 显示范围（0 ~ 0.1）
    breaks = my_breaks,         # 六段分级
    values = scales::rescale(c(0, 0.12, 0.16, 0.25))  # 4 个颜色点对应位置
  ) +
  ggtitle("M2")

ggsave('M2_Module_activity_on_UMAP.pdf',width=5, height=4.4)

edge_threshold <- 0.7  

edges <- data.frame()
for (i in 1:(nTF-1)) {
  for (j in (i+1):nTF) {
    val <- csi_mat[i,j]
    if (val > edge_threshold) {
      edges <- rbind(edges, data.frame(TF1 = tf_names[i], TF2 = tf_names[j], CSI = val))
    }
  }
}
write.csv(edges, "TF_CSI_edges_threshold0.7.csv", row.names = FALSE)
# nodes
nodes <- data.frame(TF = tf_names, Module = TF_modules)
write.csv(nodes, "TF_nodes_modules.csv", row.names = FALSE)

