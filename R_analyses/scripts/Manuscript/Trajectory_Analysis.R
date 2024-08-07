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
rm(list = ls())


load("data/RData/seurat_final.RData")
load("data/trajectory/sce_GAMed.RData")
Idents(seurat_final) <- seurat_final$celltype
DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]


Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_final))%>% unique()

nFeatures_RNA_autosomes <- colSums(seurat_final@assays$RNA$counts[Agenes,] > 0)
seurat_final <- AddMetaData(seurat_final, nFeatures_RNA_autosomes, 'nFeature_RNA_autosomes')
X_exp <- PercentageFeatureSet(seurat_final, features = Xgenes)
seurat_final <- AddMetaData(seurat_final, X_exp, 'X_exp')

Xgene_prop <- colSums(seurat_final@assays$RNA$counts[Xgenes,] > 1)/
  colSums(seurat_final@assays$RNA$counts > 1)
seurat_final <- AddMetaData(object = seurat_final, metadata = Xgene_prop, col.name = 'Xprop')


#remove outlier cells for trajectory 
DP <- subset(seurat_final, idents = c("Early cyst", "Late cyst", "Muscle"), invert = T) %>% 
  DimPlot(.)
cells1 <- CellSelector(DP, object = NULL, ident = "SelectedCells")
cells2 <- CellSelector(DP, object = NULL, ident = "SelectedCells")
cells3 <- CellSelector(DP, object = NULL, ident = "SelectedCells")
rm_cells <- c(cells1, cells2, cells3)

################################################################################
# Pseudotime analysis
sce <- sce[,!colnames(sce) %in% rm_cells]
sce_trad <- sce_trad[,!colnames(sce_trad) %in% rm_cells]

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
  mutate(`Cell type` = celltype)


cbPalette <- c("yellow3", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot() + 
  geom_point(data = Pseudotime_data, aes(x = umap_1, y = umap_2, colour = Pseudotime)) + 
  geom_path(data = traj_line[[1]], aes(x = umap_1, y = umap_2), colour = 'black', size = 1, alpha = 0.7) +
  theme_classic() +
  scale_color_continuous(low = "purple", high = "orange") + 
  labs(x = "UMAP 1", y = "UMAP 2")

p2 <- Pseudotime_data  %>% 
  ggplot(aes(x = Pseudotime, y = nFeature_RNA_autosomes, colour = `Cell type`)) + 
  geom_point(alpha = 0.5) + scale_colour_manual(values= cbPalette) +
  theme_classic() + 
  labs(y = "N\u00b0 detected autosomal genes")

p3 <- Pseudotime_data %>% 
  ggplot(., aes(x = Pseudotime, y = X_exp, colour = `Cell type`)) + 
    geom_point(alpha = 0.5) + scale_colour_manual(values= cbPalette) +
    geom_line(stat = 'smooth', aes(x =Pseudotime, y = X_exp, colour = NA),  colour = 'black', method = 'loess', se = F, 
                alpha = 0.7, size = 1) + 
    theme_classic() + 
  labs(y ="N\u00b0 expressed X genes / \n N\u00b0 expressed genes")
ggarrange(p1, p2, p3, ncol = 3)
########################## TRADESEQ ############X_exp########################## TRADESEQ ############################################

rowData(sce_trad)$assocRes <- associationTest(sce_trad, lineages = T, l2fc = log2(0.5))
assocRes <- rowData(sce_trad)$assocRes
assocRes <- assocRes[is.na(assocRes$waldStat) == F,]

assocRes %>% 
  ggplot(aes(x = pvalue_lineage1_conditionST, y = pvalue_lineage1_conditionSR)) + 
  geom_point()

bias_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionSR, "fdr") <= 0.05)
]

condRes <- conditionTest(sce_trad, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.001, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = rownames(assays(sce_trad)$counts)[oo[16]],
              alpha = 1, border = TRUE)

# least significant gene
plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = rownames(assays(sce_trad)$counts)[oo[nrow(sce_trad)]],
              alpha = 1, border = TRUE)


top_genes <- rownames(assays(sce_trad)$counts)[oo[1:20]]
filter(ortholog_table, REF_GENE_NAME %in% top_genes) %>% 
  dplyr::select(REF_GENE_NAME, consensus_gene, chr) %>% unique() %>% 
  View()

##################### MARKER GENES ---------------------------
# Look at trajectory of key marker genes for figure of in manuscript, these are 
# bb8 (gene_9679, g10101), vasa (PB.762), twe (PB.1956), and cup gene collection (STRG.2391)

# bb8
bb8 <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = "gene-9679",
              alpha = 1, border = TRUE) + 
  labs(title = "bb8", x = "")
vas <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = "PB.762",
              alpha = 1, border = TRUE) + 
  labs(title = "vasa", x = "")
twe <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = "PB.1956",
              alpha = 1, border = TRUE) + 
  labs(title = "twe")
cup <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = "STRG.2391",
              alpha = 1, border = TRUE) + 
  labs(title = "cup genes")

Fest <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = "STRG.7827",
              alpha = 1, border = TRUE) + 
  labs(title = "Fest")
CycB <- plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = "PB.4546",
              alpha = 1, border = TRUE) + 
  labs(title = "CycB")

pt_data <- data.frame(pseudotime = sce$slingPseudotime_1, celltype = sce$celltype)
pt_data$celltype <- factor(pt_data$celltype, levels = c("GSC/Spermatogonia", "Primary Spermatocytes", "Secondary Spermatocytes", "Spermatids"))
pseudotime_celltypes <- pt_data %>% 
  filter(celltype %in% c("GSC/Spermatogonia", "Primary Spermatocytes", "Secondary Spermatocytes", "Spermatids")) %>%
  ggplot(., aes(x = pseudotime, y = celltype)) + 
  geom_boxplot(outlier.alpha = 0) + 
  labs(x = "Pseudotime", y = "Cell Types") + 
  theme_bw()
# plot all together
Marker_gene_traj <- ggarrange(vas, bb8, twe, cup, Fest, CycB, p2, pseudotime_celltypes,
                              ncol = 1,
          common.legend = T, labels = c("A", "B", "C", "D", "E", "F", "G"), align = "hv")

ggsave("plots/Marker_gene_traj.pdf", Marker_gene_traj, width = 10, height = 15)


##################### look into ###############
yhatSmooth <- predictSmooth(sce_trad, gene = c("gene-9128", "PB.370"), nPoints = 50, tidy = FALSE)
yhatSmooth <- yhatSmooth[is.na(yhatSmooth[,1]) == FALSE,]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = T, 
                       show_colnames = T)


####################################################



########### MAYBE CHANGE? -----------


# test for dynamic expression of X chromosome ---

counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)
counts_exp <- rbind(counts[gois$geneID[gois$celltype == "Spermatids"],], 
                    seurat_marker@meta.data$xprop)

X_exp <- counts(sce[Xgenes,]) %>% colSums()

rownames(counts_exp)[nrow(counts_exp)] <- "X expression"

rm_cells <- which(cellWeights ==0)

sce_trad <- fitGAM(#counts = counts[gois$geneID[gois$celltype == "Spermatids"],], 
  counts = counts_exp[,-rm_cells],
  pseudotime = pseudotime[-rm_cells,1], 
  cellWeights = cellWeights[-rm_cells,1],
  conditions = factor(colData(sce)$treatment[-rm_cells]),
  nknots = 5, parallel=F)


rowData(sce_trad)$assocRes <- associationTest(sce_trad, lineages = T, l2fc = log2(0.5))
assocRes <- rowData(sce_trad)$assocRes
assocRes <- assocRes[is.na(assocRes$waldStat) == F,]

assocRes %>% 
  ggplot(aes(x = pvalue_lineage1_conditionST, y = pvalue_lineage1_conditionSR)) + 
  geom_point()

bias_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionSR, "fdr") <= 0.05)
]

yhatSmooth <- predictSmooth(sce_trad, gene = c("X expression"), nPoints = 50, tidy = FALSE)
yhatSmooth <- yhatSmooth[is.na(yhatSmooth[,1]) == FALSE,]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = T, 
                       show_colnames = T)
