################# load libraries and functions and data prep ###################
library(dplyr)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(ggrepel)
library(lme4)
library(GenomicFeatures)
library(ggrepel)
library(lme4)
rm(list = ls())

load("data/RData/seurat_final.RData")
DefaultAssay(seurat_final) <- "RNA"
seurat_final$treatment <- "SR"
seurat_final$treatment[grep("st", seurat_final$sample)] <- "ST"
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene, FBconcensus) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

seurat_final <- JoinLayers(seurat_final)

################################################################################


######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
n=1000
replicates <- list()
counts <- seurat_final@assays$RNA$counts %>% as.data.frame()
metadata <- seurat_final@meta.data
ct = "Muscle"
write.table(metadata, file = 'data/pseudoreplicate/metadata_pr.csv', quote = F, row.names = T, sep = ',')
write.table(counts, file = 'data/pseudoreplicate/counts_pr.csv', quote = F, row.names = T, sep = ',')
i=1
for (i in 1:n){
  print(i)
  dataset <- list()
  for (ct in unique(metadata$celltype)){
    
    print(ct)
    cl_c=0
    cells_all <- list()
    for (trt in c("ST", "SR")){
      cl_c=cl_c+1
      cells <- c()
      for (s in 1:4){
        ps <- sample(metadata$cells[metadata$celltype == ct & metadata$treatment == trt], 100, replace = T)
        cells <-unique(c(cells,ps)) 
        
        s1 <- counts[,ps] %>% rowSums()
        if (s == 1 & trt == "ST"){
          s1_counts <- s1
        } else {
         s1_counts <- cbind(s1_counts, s1)
        }
        
      }
      cells_all[[cl_c]] <- cells
      
      
    }
    g1_keep <- (rowSums(counts[,cells_all[[1]]] > 1) > (0.05 * length(cells_all[[1]]))) #| 
    #rowSums(seurat_final@assays$RNA$counts[,ST_cells] > 1) > 10)
    
    g2_keep <- (rowSums(counts[,cells_all[[2]]] > 1) > (0.05 * length(cells_all[[2]])))# | 
    #                  rowSums(seurat_final@assays$RNA$counts[,SR_cells] > 1) > 10)
    
    keep <- g1_keep | g2_keep
    colnames(s1_counts) <- c(paste0("ST", 1:4), paste0("SR", 1:4))
    dataset[[ct]] <- s1_counts[keep,]
  }
  replicates[[paste0("rep", i)]] <- dataset
}


c=0
results_df <- NULL
for (i in 1:n){
  print(i)
  for (ct in unique(metadata$celltype)){
    print(ct)
    c=c+1
    counts <- replicates[[paste0("rep", i)]][[ct]]
    
    
    y <- DGEList(counts, samples = colnames(counts), group = c(rep("ST", 4), rep("SR", 4)))
    y <- calcNormFactors(y)
    design <- model.matrix(~group, y$samples)
    y <- normLibSizes(y)
    y <- estimateDisp(y, design)
    cpm <- edgeR::cpm(y, log = F)
    
    keep1 <- cpm[,y$samples$group== "ST"] >= 1
    keep1 <- rowSums(keep1) > (ncol(keep1)/2)
    keep2 <- cpm[,y$samples$group== "SR"] >= 1
    keep2 <- rowSums(keep2) > (ncol(keep2)/2)
    keep <- keep1 | keep2
    
    y <- y[keep,keep.lib.sizes=FALSE]
      
    fit <- glmQLFit(y, design, robust = T)
    res <- glmQLFTest(fit, coef=ncol(design))
    toptags <- topTags(res, n = nrow(res))
    out <- toptags$table %>% 
      mutate(gene = rownames(.), 
             celltype = ct, 
             rep = paste0("rep", i), 
             treatment = "ST")
    if (c == 1){
      results_df <- out
    } else {
      results_df <- rbind(results_df, out)
    }
  }
}



write.table(results_df[1:100,], file = "del.csv", sep = ",", row.names = F, quote = F)

STST

results_df %>% 
  mutate(sig = ifelse(FDR < 0.05 & abs(logFC) >=1, "Sig", "Not")) %>% 
  group_by(celltype, rep, sig) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = celltype, y = n, fill = sig)) + geom_boxplot()









