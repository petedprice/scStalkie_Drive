#### X_LINKED CELLS ----
library(GenomicFeatures)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(ggpubr)
load("data/RData/marker_seurat.RData")

#Tidying up marker data 
seurat_marker$celltype <- seurat_marker$sctype_labels
seurat_marker$celltype <- gsub("_", " ", seurat_marker$celltype)
seurat_marker$celltype <- gsub(" or ", ", ", seurat_marker$celltype)
seurat_marker$treatment <- "ST"
seurat_marker$treatment[grepl("^sr", seurat_marker$sample)] <- "SR"

#repeat for seurat_marker_noX
seurat_marker_noX$celltype <- seurat_marker_noX$sctype_labels
seurat_marker_noX$celltype <- gsub("_", " ", seurat_marker_noX$celltype)
seurat_marker_noX$celltype <- gsub(" or ", ", ", seurat_marker_noX$celltype)
seurat_marker_noX$treatment <- "ST"
seurat_marker_noX$treatment[grepl("^sr", seurat_marker_noX$sample)] <- "SR"


ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

Xgenes <- filter(ortholog_table, chr == "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  c()


Agenes <- filter(ortholog_table, chr != "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>%
  c()
propX <- length(Xgenes)/length(genes)

#Plots
plots <- list()
plots[[1]] <- DimPlot(seurat_marker, group.by = 'celltype', split.by = 'treatment', label = T) + 
  labs(title = "X linked genes included")
plots[[2]] <- DimPlot(seurat_marker_noX, group.by = 'celltype', split.by = 'treatment', label = T) + 
  labs(title = "X linked genes removed")
arranged <- ggarrange(plotlist = plots) 

ggsave("plots/X_linked_exp_comp.pdf", arranged, height = 10, width = 25)


###### CALCULATING X LINKED EXPRESSION AND ADDING TO OBJECT ---------
counts <- GetAssayData(seurat_marker, assay = "RNA")
Total_expression <- colSums(counts)
X_linked_expression <- counts[
  rownames(counts) %in% Xgenes,] %>% 
  colSums()
X_prop <- (X_linked_expression/Total_expression)

seurat_marker <- AddMetaData(seurat_marker, X_prop, 'xprop')
Ameans <- colMeans(counts[Agenes[Agenes %in% rownames(counts)],])
Xmeans <- colMeans(counts[Xgenes[Xgenes %in% rownames(counts)],])
means <- colMeans(counts)
ratios <- Ameans/Xmeans
#add data to metadata 
seurat_marker <- AddMetaData(seurat_marker, ratios, 'ratios')
#add Ameans and Xmeans
seurat_marker <- AddMetaData(seurat_marker, Ameans, 'Ameans')
seurat_marker <- AddMetaData(seurat_marker, Xmeans, 'Xmeans')
seurat_marker <- AddMetaData(seurat_marker, means, 'means')

#seurat_marker$celltype[seurat_marker$seurat_clusters == 22] <- "Y-sperm"


plots <- list()
#plots[[1]] <- 
seurat_marker@meta.data %>% 
  ggplot(aes(x = treatment, y = xprop, fill = celltype, colour = sample)) +
  facet_grid(~celltype) +
  geom_boxplot() 


results <- list()
for (c in unique(seurat_marker$celltype)){
  results[[c]] <- seurat_marker@meta.data %>% 
    filter(celltype == c) %>%
    lmerTest::lmer(xprop ~ treatment + (1|sample), .) %>% 
    summary()
}

plots[[2]] <- DimPlot(seurat_marker, group.by = 'seurat_clusters', label = T, split.by = 'treatment')
 
plots[[3]] <- ggplot(seurat_marker@meta.data, aes(x = xprop, fill = celltype)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~treatment) + 
  labs(title = "X-linked expression as a proportion of total expression by celltype")

plots[[4]] <- ggplot(seurat_marker@meta.data, aes(x = celltype, y = xprop, fill = celltype)) + 
  geom_boxplot() + 
  facet_wrap(~treatment) + 
  labs(title = "X-linked expression as a proportion of total expression by celltype")

arranged <- ggarrange(plotlist = plots, nrow = 4)
ggsave("plots/findingY.pdf", arranged, height = 20, width = 20)


save(seurat_marker, file = "data/RData/findingY.RData")


