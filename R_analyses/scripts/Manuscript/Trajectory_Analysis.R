########################## DATA AND PACKAGES ################################### 
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
load("data/RData/DEG_DC.RData")

Idents(seurat_final) <- seurat_final$celltype
DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

nFeatures_RNA_autosomes <- colSums(seurat_final@assays$RNA$counts[Agenes,] > 1)
seurat_final <- AddMetaData(seurat_final, nFeatures_RNA_autosomes, 'nFeature_RNA_autosomes')
X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  colSums(seurat_final@assays$RNA$counts > 1)
seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')

cells2 <- colnames(sce)
seurat_final <- seurat_final %>% 
  subset(cells = cells2)


################################################################################
# Pseudotime analysis
Pseudotime_data <- reducedDims(sce)$UMAP %>% as.data.frame() %>% 
  merge(seurat_final@meta.data, by.x = 0, by.y = 'cells', all.y = F) %>% 
  mutate(cells = Row.names)

sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC/Spermatogonia',
                 end.clus = 'Late spermatids')


traj_line <- SlingshotDataSet(sce)@curves$Lineage1

pts <- data.frame(Pseudotime = c(unlist(traj_line[[3]])), 
                  cells = names(c(unlist(traj_line[[3]]))))

Pseudotime_data <- Pseudotime_data %>% 
  merge(pts, by = 'cells') %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate('Cell type' = celltype)



########################## TRADESEQ ############################################

#Gene filter #
ST_cells <- filter(seurat_final@meta.data, treatment == "ST")$cells %>% intersect(colnames(sce_trad))
SR_cells <- filter(seurat_final@meta.data, treatment == "SR")$cells %>% intersect(colnames(sce_trad))

ST_keep <- rowSums(counts(sce_trad[,ST_cells]) > 1) > (0.1 * length(ST_cells))
SR_keep <- rowSums(counts(sce_trad[,SR_cells]) > 1) > (0.1 * length(SR_cells))
keep <- ST_keep | SR_keep

sce_trad <- sce_trad[keep,]
rowData(sce_trad)$assocRes <- associationTest(sce_trad, lineages = T, l2fc = log(2))
assocRes <- rowData(sce_trad)$assocRes

assocRes$ST_padj <- p.adjust(assocRes$pvalue_lineage1_conditionST, "fdr")
assocRes$SR_padj <- p.adjust(assocRes$pvalue_lineage1_conditionSR, "fdr")

assocRes <- filter(assocRes, ST_padj < 0.05 | SR_padj < 0.05)
assocRes <- assocRes[is.na(assocRes$waldStat) == F,]

sce_trad <- sce_trad[rownames(assocRes),]

gois1 <- rownames(assocRes)[assocRes$ST_padj <=0.05] %>% .[!is.na(.)]

smooth_ST_genes <- predictSmooth(sce_trad, gene = gois1, nPoints = 100, tidy = F)
heatmap_ST <- pheatmap::pheatmap(t(scale(t(smooth_ST_genes[, 1:50]))),
                                      cluster_cols = FALSE,
                                      show_rownames = FALSE, 
                                      show_colnames = FALSE)


condRes <- conditionTest(sce_trad, l2fc = 2) %>% 
  mutate(gene = rownames(.)) %>% 
  merge(ortholog_table, by.x = 'gene', by.y = 'REF_GENE_NAME', all.y = F)

condRes$padj <- p.adjust(condRes$pvalue, "fdr")
condRes <- filter(condRes, padj < 0.05 & !is.na(waldStat)) 

oo <- order(condRes$waldStat, decreasing = TRUE)
condRes %>% 
  mutate(consensus_gene = ifelse(consensus_gene == gene, " ", consensus_gene)) %>%
  rename(`Teleopsis dalmanni Gene` = gene, 
         `Chr` = chr, `Drosophila ortholog` = consensus_gene, `FDR-adjusted p-value` = padj
         ) %>% 
  .[,c("Teleopsis dalmanni Gene", "Drosophila ortholog", "Chr", "waldStat", "df", 
       "pvalue", "FDR-adjusted p-value")] %>% 
  
  write.table(., "data/MANUSCRIPT/differential_trajectories.tsv", sep = "\t", quote = F, row.names = F)

# most significant gene
Smoother_func <- function(sce_trad, gn, ot = ortholog_table){
  chr = filter(ot, REF_GENE_NAME == gn)$chr
  if(chr == "Chr_X"){chr = " (X)"}else{chr = " (Autosome)"}
  title <- filter(ot, REF_GENE_NAME == gn)%>% 
    dplyr::select(consensus_gene, REF_GENE_NAME) %>% unique()
  if(title$REF_GENE_NAME != title$consensus_gene){
    title = paste0("D.mel: ", title$consensus_gene, chr)
  } else {
    title = paste0("Ref gene: ", title$REF_GENE_NAME, chr)
  }
  plot <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
                        gene = gn,
                        alpha = 1, border = TRUE, curvesCols = c("#FC8D62", "#66C2A5")
                        ) + 
    ggtitle(title) +
    scale_color_manual(values = c("#FC8D62", "#66C2A5"), labels = c("SR", "ST"), name = "") + 
    guides(colour = guide_legend(override.aes = list(size = 2))) + 
    theme(legend.position = 'bottom', 
          legend.text = element_text(size = 15)) + labs(x = "", y = "")
  
}

# GO term enrichment 

background <- filter(ortholog_table, REF_GENE_NAME %in% rownames(assocRes)) 
GO_data <- condRes 

ego <- enrichGO(gene          = unique(GO_data$consensus_gene), 
                  universe = unique(background$consensus_gene),
                OrgDb         = org.Dm.eg.db, 
                keyType = 'SYMBOL',
                ont           = "ALL",
                pAdjustMethod = "bonferroni",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05, readable= TRUE)


ego <- ego %>% #round all numeric value to 3dp
  dplyr::mutate(across(where(is.numeric), round, 4))


write.table(ego, "data/MANUSCRIPT/TRADESEQ_GO_terms.tsv", sep = "\t", quote = F, row.names = F)
system("open data/MANUSCRIPT/TRADESEQ_GO_terms.tsv")

############## CANDIDATE GENES ############

DGE_DI <- merge(dif_exp_data, condRes, by.x = 'genes', by.y = 'gene', all = F) %>% 
  filter(Significant != 'Unbiased') %>% 
  filter(grepl("perm", celltype)) %>% 
  dplyr::select(waldStat, genes, consensus_gene.x) %>% unique() %>% 
  #filter(genes != consensus_gene.x) %>% 
  top_n(20, waldStat) %>% 
  arrange(desc(waldStat))

gois <- DGE_DI$genes %>% unique()

plots <- lapply(gois, Smoother_func, sce_trad = sce_trad)
arranged <- ggarrange(plotlist = plots, ncol =4, nrow = 5, common.legend = T) + 
  theme(legend.text = element_text(size = 30), # Change legend text size
        legend.title = element_text(size = 30))
arranged <- annotate_figure(arranged, bottom = text_grob("Pseudotime", size = 14, vjust = -0.5), left = text_grob("log(expression + 1)", , size = 14, rot = 90, vjust = 1))
ggsave("plots/S9_TRADESEQ2_top20_genes.pdf", arranged, width = 15, height = 12)
system("open plots/S9_TRADESEQ2_top20_genes.pdf")


################################################################################
################################################################################

save(Pseudotime_data, Smoother_func, traj_line, file = 'data/RData/traj_analyis.RData')

################################################################################
################################################################################




























# 
# 
# 
# ### CELL TYPE ABUNDANCE AGAINST TIME 
# 
# p1 <- Pseudotime_data %>% 
#   ggplot(aes(x = Pseudotime, fill = treatment)) + 
#   geom_density(alpha = 0.5) +scale_color_manual(values = c("#FC8D62", "#66C2A5"), labels = c("SR", "ST"), name = "") + 
#   theme_classic()
# out <- cowplot::plot_grid(p1, p4, nrow = 2, align = 'v', axis = 'lr')
# ggsave("plots/TRADESEQ2_pseudotime_density.pdf", out, width = 6, height = 3)
# system("open plots/TRADESEQ2_pseudotime_density.pdf")
# ##############################################################################
# ##########################WORK IN PROGRESS####################################
# ##############################################################################
# 
# condRestemp <- condRes #%>% filter(chr == "Chr_X")
# yhatSmooth <- predictSmooth(sce_trad, gene = condRestemp$gene, nPoints = 50, tidy = FALSE)
# 
# yhatSmoothScaled <- t(scale(t(yhatSmooth)))
# 
# color_palette <- colorRampPalette(c("green", "white", "purple"))(50)
# breaks <- seq(-6, 6, length.out = 51)
# 
# 
# heatSmooth_ST <- pheatmap::pheatmap(yhatSmoothScaled[, 51:100],
#                            cluster_cols = FALSE,
#                            show_rownames = FALSE, show_colnames = FALSE, main = "ST", legend = T, 
#                            silent = T#, #color = color_palette,      # Use the custom color palette
#                           #breaks = breaks,            # Apply the custom breaks
#                           #legend_breaks = seq(-6, 6, 1),  # Set the legend breaks for clarity
#                           #legend_labels = seq(-6, 6, 1)
# )
# 
# order <- heatSmooth_ST$tree_row$order
# heatSmooth_SR <- pheatmap::pheatmap(yhatSmoothScaled[heatSmooth_ST$tree_row$order, 1:50],
#                                  cluster_cols = FALSE, cluster_rows = FALSE,
#                                  show_rownames = FALSE, show_colnames = FALSE, main = "SR",
#                                  legend = T, silent = T #, #color = color_palette,      # Use the custom color palette
# #                           breaks = breaks,            # Apply the custom breaks
# #                           legend_breaks = seq(-6, 6, 1),  # Set the legend breaks for clarity
# #                           legend_labels = seq(-6, 6, 1)
# )
# 
# heatSmooth_SRST <- pheatmap::pheatmap(yhatSmoothScaled[, 1:50][,]-yhatSmoothScaled[, 51:100][,],
#                           cluster_cols = FALSE, cluster_rows = T,
#                           show_rownames = FALSE, show_colnames = FALSE, main = "SR-ST",
#                           legend = TRUE, silent = TRUE#, 
#                          # color = color_palette,      # Use the custom color palette
#                          # breaks = breaks,            # Apply the custom breaks
#                         #  legend_breaks = seq(-2.5, 6, 1),  # Set the legend breaks for clarity
#                          # legend_labels = seq(-2.5, 6, 1)
# )
# 
# dev.off()
# 
# gridExtra::grid.arrange(heatSmooth_ST[[4]], heatSmooth_SR[[4]], heatSmooth_SRST[[4]],ncol = 3)
# a <- ggplotGrob(heatSmooth_ST)
# pt_ct <- Pseudotime_data %>% 
#   ggplot(aes(x = Pseudotime, y = celltype)) + geom_boxplot()
# dev.off()
# cowplot::plot_grid(heatSmooth_ST, pt_ct, nrow = 2, align = 'hv', axis = 'rltb')
# grid.arrange(heatSmooth_ST, pt_ct, nrow = 2)
# ##################### MARKER GENES ---------------------------
# 
# # Look at trajectory of key marker genes for figure of in manuscript, these are 
# # bb8 (gene_9679, g10101), vasa (PB.762), twe (PB.1956), and cup gene collection (STRG.2391)
# 
# # bb8
# bb8 <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
#               gene = "gene-9679",
#               alpha = 1, border = TRUE) + 
#   labs(title = "bb8", x = "")
# vas <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
#               gene = "PB.762",
#               alpha = 1, border = TRUE) + 
#   labs(title = "vasa", x = "")
# twe <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
#               gene = "PB.1956",
#               alpha = 1, border = TRUE) + 
#   labs(title = "twe")
# cup <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
#               gene = "STRG.2391",
#               alpha = 1, border = TRUE) + 
#   labs(title = "cup genes")
# 
# Fest <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
#               gene = "STRG.7827",
#               alpha = 1, border = TRUE) + 
#   labs(title = "Fest")
# CycB <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
#               gene = "PB.4546",
#               alpha = 1, border = TRUE) + 
#   labs(title = "CycB")
# 
# pt_data <- data.frame(pseudotime = sce$slingPseudotime_1, celltype = sce$celltype)
# pt_data$celltype <- factor(pt_data$celltype, levels = c("GSC/Spermatogonia", "Primary Spermatocytes", "Secondary Spermatocytes", "Spermatids"))
# pseudotime_celltypes <- pt_data %>% 
#   filter(celltype %in% c("GSC/Spermatogonia", "Primary Spermatocytes", "Secondary Spermatocytes", "Spermatids")) %>%
#   ggplot(., aes(x = pseudotime, y = celltype)) + 
#   geom_boxplot(outlier.alpha = 0) + 
#   labs(x = "Pseudotime", y = "Cell Types") + 
#   theme_bw()
# # plot all together
# Marker_gene_traj <- ggarrange(vas, bb8, cup, Fest, CycB, p2, pseudotime_celltypes,
#                               ncol = 1,
#           common.legend = T, labels = c("A", "B", "C", "D", "E", "F", "G"), align = "hv")
# 
# ggsave("plots/Marker_gene_traj.pdf", Marker_gene_traj, width = 10, height = 15)
# system("open plots/Marker_gene_traj.pdf")
# ####################################################`