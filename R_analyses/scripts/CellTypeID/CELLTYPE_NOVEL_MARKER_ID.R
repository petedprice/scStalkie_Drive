library(Seurat)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(slingshot)
library(tradeSeq)
library(grDevices)
library(TSCAN)
library(pheatmap)
library(gridExtra)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Dm.eg.db)
library(cowplot)


rm(list = ls())


load("data/RData/seurat_final.RData")
load("data/trajectory/sce_GAMed.RData")
load("data/RData/DEG_DC.RData")``

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))


Idents(seurat_final) <- "celltype"
germ <- seurat_final$celltype %>% unique() %>% .[grepl("perm",.)]
seurat_final <- subset(seurat_final, idents = germ) 
Idents(seurat_final) <- "integrated_snn_res.0.4"

### FINDALLMARKERS APPROACH ###

out <- list()
for (ct in unique(seurat_final$celltype)){
  markers <- seurat_final %>% 
    FindMarkers(only.pos = T, ident.1 = ct) %>% 
    filter(p_val_adj == 0 & avg_log2FC > 1) %>% 
    mutate(gene = row.names(.))
  
  seurat_subset <- subset(seurat_final, idents = ct)
  keep_genes <- rowSums(seurat_subset@assays$RNA@counts > 1) > (0.05 * ncol(seurat_subset))
  background <- names(keep_genes)[keep_genes]
  sskm <- markers$gene
  markers_dros = filter(ortholog_table, REF_GENE_NAME %in% sskm)$consensus_gene %>% unique()
  background_dros = filter(ortholog_table, REF_GENE_NAME %in% background)$consensus_gene %>% unique()
  
  
  ego <- enrichGO(gene          = markers_dros, 
                  universe = background_dros,
                  OrgDb         = org.Dm.eg.db, 
                  keyType = 'SYMBOL',
                  ont           = "ALL",
                  pAdjustMethod = "bonferroni",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05, readable= TRUE)
  
  ego <- as.data.frame(ego) %>% 
    mutate(cluster = ct)
  
  out[[ct]] = ego
}

enrich_data <- do.call(rbind, out)
View(enrich_data)
