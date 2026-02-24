rm(list = ls())
library(psych)
library(tidyverse)
library(Seurat)
library(dplyr)
library(openxlsx)
library(ggplot2)

setwd('D:/Desktop')
dat = readRDS('mnn.rds')

dat = dat[,dat@meta.data$celltype %in% c('Astro', 'BG', 'Epen','Micro','Neu','OEC')]
dat = dat[,dat@meta.data$celltype %in% c('OPC', 'PC', 'Pit', 'PVM', 'TC', 'VEC','VLMC','VSMC')]
dat = dat[,dat@meta.data$celltype %in% c('Oligo')]

dat@meta.data$celltype[dat@meta.data$celltype == 'BG'] <- 'Astro'

mat <- GetAssayData(dat, layer = "data")
mat <- as.data.frame(t(as.matrix(mat)))
info <- dat@meta.data
info2 <- as.data.frame(info[,c(16,29,30)])
rownames(info2)=rownames(info)
dat <- merge(mat,info2,by="row.names",sort=FALSE)
rownames(dat)<-dat[,1]
dat<- dat[-1]
dim(dat)

dat2=dat
dat2$Organ_celltype = paste(dat2$Organ,dat2$celltype,sep = "_")

table_counts <- table(dat2$Organ_celltype)

keep_groups <- names(table_counts[table_counts >= 100])

dat2 <- subset(dat2, subset = Organ_celltype %in% keep_groups)

celltypes <- unique(dat2$celltype)
pvalue1 <- NULL

for (ct in celltypes) {
  message("Processing celltype: ", ct)
  
  all_results <- data.frame()
  all_lm_stats <- list()  # 用于保存每个细胞类型 + 器官下的完整 t/p 值矩阵
  
  dat_ct <- dat %>% filter(celltype == ct)

  unique_sexes <- unique(na.omit(dat_ct$Sex1))
  num_sex <- length(unique_sexes)
  
  if (num_sex > 1) {
    dat_ct$Sex1 <- factor(dat_ct$Sex1)
    model_formula_str <- "expr ~ Group + Sex1"
    message(" - [Info] Multiple sexes detected. Adjusting for Sex.")
  } else {
    
    model_formula_str <- "expr ~ Group"
    message(" - [Info] Single sex detected. Removing Sex from model.")
  }
  
  dat_ct$Organ <- factor(dat_ct$Organ)
  
  
  if (length(unique(dat_ct$Organ)) == 1) {
    message(" - Skipping celltype ", ct, " (only one organ available)")
    next  
  }
  
  organs <- c('Brain stem',
              'Cerebellum',
              'Cerebral cortex',
              'Cingulate gyrus',
              'Corpus callosum',
              'Hippocampus',
              'Hypothalamus',
              'Olfactory bulb',
              'Spinal cord',
              'Striatum',
              'Thalamus')
  
  for (target_organ in organs) {
    if (!(target_organ %in% dat_ct$Organ)) {
      message(" - Skipping organ (no data): ", target_organ)
      next
    }
    message(" - Analyzing organ: ", target_organ)
    
    dat_ct$Group <- ifelse(dat_ct$Organ == target_organ, target_organ, "Other")
    dat_ct$Group <- factor(dat_ct$Group, levels = c("Other", target_organ))
    
    pvals <- c()
    estimates <- c()
    upregulated <- 0
    downregulated <- 0
    pvalue1 <- data.frame()  
    
    
    for (i in 1:27849) {
      
      expr <- dat_ct[[i]]
      fit <- try(lm(as.formula(model_formula_str), data = dat_ct), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        coef_table <- summary(fit)$coefficients
        
        row_target <- ifelse("Group" %in% rownames(coef_table),
                             "Group",
                             paste0("Group", target_organ))
        
        if (row_target %in% rownames(coef_table)) {
          est <- coef_table[row_target, "Estimate"]
          p <- coef_table[row_target, "Pr(>|t|)"]
          
          estimates <- append(estimates, est)
          pvals <- append(pvals, p)
          
          if (est > 0 && p < 0.01 / 27849) {
            upregulated <- upregulated + 1
          } else if (est < 0 && p < 0.01 / 27849) {
            downregulated <- downregulated + 1
          }
        }
        
        result <- as.data.frame(coef_table)[, c("Estimate", "t value", "Pr(>|t|)")]
        result <- as.data.frame(t(result))
        
        est_res <- result[1, ]
        colnames(est_res) <- paste0("Estimate_", colnames(est_res))
        
        tres <- result[2, ]
        colnames(tres) <- paste0("T_", colnames(tres))
        
        pres <- result[3, ]
        colnames(pres) <- paste0("P_", colnames(pres))
        
        tmp1 <- cbind(est_res, tres, pres)
        rownames(tmp1) <- colnames(dat_ct)[i]  
        
        pvalue1 <- rbind(pvalue1, tmp1)
      }
    }
    
    deg_count <- sum(pvals < 0.01 / 27849, na.rm = TRUE)
    

    all_results <- rbind(all_results,
                         data.frame(Celltype = ct,
                                    Organ = target_organ,
                                    DEG_count = deg_count,
                                    Upregulated = upregulated,
                                    Downregulated = downregulated))
    
   
    all_lm_stats[[paste(ct, target_organ, sep = "_")]] <- pvalue1
  }
  write.csv(all_results, 
            file = paste0("summary_", ct, ".csv"), 
            row.names = FALSE)
  
  ct_lm_stats <- all_lm_stats[grepl(paste0("^", ct, "_"), names(all_lm_stats))]
  
  for (name in names(ct_lm_stats)) {
    
    safe_name <- make.names(name)
    
    file_path <- paste0("lm_stats_", safe_name, ".csv")
    
    write.csv(ct_lm_stats[[name]], 
              file = file_path, 
              row.names = TRUE)
  }
}


