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
Idents(seurat_final) <- seurat_final$celltype

keep_clusters <- c("GSC/Spermatogonia", "Primary spermatocytes", "Spermatocytes", "Secondary spermatocytes", "Spermatids")
keep_clusters <- keep_clusters[keep_clusters %in% unique(seurat_final$celltype)]

Idents(seurat_final) <- 'integrated_snn_res.0.4'

sce <- seurat_final %>% 
  subset(., idents = keep_clusters) %>% 
  as.SingleCellExperiment(., assay = 'RNA')

sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC/Spermatogonia',
                 end.clus = 'Spermatids')

counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)


ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

################################################################################

Pseudotime_data <- reducedDims(sce)$UMAP %>% as.data.frame() %>% 
  merge(seurat_final@meta.data, by.x = 0, by.y = 'cells', all.y = F) %>% 
  mutate(cells = Row.names)


sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC/Spermatogonia',
                 end.clus = 'Spermatids')


traj_line <- SlingshotDataSet(sce)@curves$Lineage1

pts <- data.frame(Pseudotime = c(unlist(traj_line[[3]])), 
                  cells = names(c(unlist(traj_line[[3]]))))

Pseudotime_data <- Pseudotime_data %>% 
  merge(pts, by = 'cells')

p1 <- ggplot() + 
  geom_point(data = Pseudotime_data, aes(x = umap_1, y = umap_2, colour = Pseudotime)) + 
  geom_smooth(data = traj_line[[1]], aes(x = umap_1, y = umap_2), colour = 'azure4') +
  theme_bw() +
  scale_color_continuous(low = "purple", high = "orange") + 
  labs(x = "UMAP 1", y = "UMAP 2")

Pseudotime_data  %>% 
  ggplot(aes(x = Pseudotime, y = nFeature_RNA)) + 
  geom_point() + 
  geom_smooth()
seurat_final@meta.data %>% 
  ggplot(aes(x = nFeature_RNA, y = log10GenesPerUMI, colour = Phase)) + 
  geom_point()
########################## TRADESEQ ############################################

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
              gene = rownames(assays(sce_trad)$counts)[oo[17]],
              alpha = 1, border = TRUE)

# least significant gene
plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = rownames(assays(sce_trad)$counts)[oo[nrow(sce_trad)]],
              alpha = 1, border = TRUE)


top_genes <- rownames(assays(sce_trad)$counts)[oo[1:20]]
View(filter(ortholog_table, REF_GENE_NAME %in% top_genes))

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
Marker_gene_traj <- ggarrange(vas, bb8, twe, cup, Fest, CycB, pseudotime_celltypes,
                              ncol = 1,
          common.legend = T, labels = c("A", "B", "C", "D", "E", "F"), align = "hv")

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
Xgenes <- filter(ortholog_table, chr == "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  c() %>% intersect(rownames(seurat_integrated_ss))

DefaultAssay(seurat_integrated_ss) <- "RNA"
Xexp <- PercentageFeatureSet(seurat_integrated_ss, features = Xgenes)
seurat_integrated_ss <- AddMetaData(object = seurat_integrated_ss, metadata = Xexp, col.name = 'Xexp')

Xgene_prop <- colSums(seurat_integrated_ss@assays$RNA$counts[Xgenes,] > 0)/
  colSums(seurat_integrated_ss@assays$RNA$counts > 0)

seurat_integrated_ss <- AddMetaData(object = seurat_integrated_ss, metadata = Xgene_prop, col.name = 'Xprop')

#######################################


#######################
Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype
sce <- seurat_integrated_ss %>% 
  subset(., idents = c("GSC & Spermatogonia", "Spermatocytes", "Spermatids")) %>% 
  as.SingleCellExperiment(., assay = 'RNA')

sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC & Spermatogonia',
                 end.clus = 'Spermatids')

counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)
colorspt <- colorRampPalette(brewer.pal(11,'RdBu')[-6])(100)
plotcolpt <- colorspt[cut(sce$slingPseudotime_1, breaks=100)]

colorsct <- brewer.pal(6,'Set2')
names(colorsct) <- unique(sce$celltype)
colorsct2 <- brewer.pal(11,'Set3')
names(colorsct2) <- unique(sce$integrated_snn_res.0.1)
#use a brewer.pall pallete with 16 colours for colorsct2



plotcolct <- colorsct[sce$celltype]
plotcolct2 <- colorsct2[sce$integrated_snn_res.0.1]
dev.off()
layout(matrix(1:3, nrow = 1))
par(mar = c(4.5,4,1,1))
plot(reducedDims(sce)$UMAP, col = plotcolct, pch=16, asp = 1)
plot(reducedDims(sce)$UMAP, col = plotcolct2, pch=16, asp = 1)
plot(reducedDims(sce)$UMAP, col = plotcolpt, pch=16, asp = 1)
lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black')

seurat_integrated_ss <- AddMetaData(seurat_integrated_ss, pseudotime[,1], col.name = 'pseudotime')
plots <- list()
plots[[1]] <- seurat_integrated_ss@meta.data %>% 
  ggplot(aes(x = pseudotime, y = Xprop, colour = treatment)) + #geom_point(alpha = 0.1) + 
  geom_smooth() + 
  theme_classic() + 
  scale_colour_brewer('Greens')

plots[[2]] <- seurat_integrated_ss@meta.data %>% 
  ggplot(aes(x = pseudotime, colour = treatment)) + geom_density() + 
  scale_colour_brewer('Greens') + 
  theme_classic()

plots[[3]] <- ggplot(seurat_integrated_ss@meta.data, aes(x = pseudotime, y = celltype)) + 
  geom_boxplot(outlier.alpha = 0)

ggarrange(plotlist = plots, nrow = 3, common.legend = TRUE, legend = "bottom")




figure <- ggarrange(plots[[1]]  + rremove("xlab"), plots[[2]] + rremove("xlab"), plots[[3]] + rremove("xlab"),# remove axis labels from plots
                    labels = NULL,
                    nrow = 3,
                    common.legend = TRUE, legend = "bottom",
                    align = "hv")


figure
seurat_integrated_ss@meta.data %>% 
  filter(celltype %in% c("GSC & Spermatogonia", "Spermatocytes", "Spermatids")) %>%
  ggplot(., aes(x = pseudotime, y = Xprop, color = treatment)) +
    geom_point(size = 0.5) +
    #geom_smooth(method = 'loess', se = FALSE, span = 0.1) +
    scale_color_manual(values = c('red', 'blue')) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(x = 'Pseudotime', y = 'X-linked gene expression')
#    facet_wrap(~celltype, ncol = 4)



################################# TRADESEQ -------------------------------------

gois <- lapply(cell_type_DEG_output, function(x)(return(rbind(x$ST_bias, 
                                                              x$SR_bias)))) %>% 
  bind_rows()

#evaluate k 
icMat <- evaluateK(sce,
                   conditions = factor(sce$treatment),
                   nGenes = 200,
                   k = 3:10, parallel=F)

#FItting GAM
sce <- fitGAM(counts = counts(sce), pseudotime = pseudotime, 
              cellWeights = cellWeights, nknots = 4, conditions = as.factor(sce$treatment))





 
#################################################





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


condRes <- conditionTest(sce_trad, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = rownames(assays(sce_trad)$counts)[oo[1]],
              alpha = 1, border = TRUE)

# least significant gene
plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = rownames(assays(sce_trad)$counts)[oo[nrow(sce_trad)]],
              alpha = 1, border = TRUE)
