library(Seurat)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

load("data/RData/marker_seurat.RData")

#Tidying up marker data 
seurat_marker$celltype <- seurat_marker$sctype_labels
seurat_marker$celltype <- gsub("_", " ", seurat_marker$celltype)
seurat_marker$celltype <- gsub(" or ", ", ", seurat_marker$celltype)
seurat_marker$treatment <- "ST"
seurat_marker$treatment[grepl("^sr", seurat_marker$sample)] <- "SR"

pg_markers <- readxl::read_excel("data/PopGroup_markers.xlsx") 
colnames(pg_markers)[1] <- "seurat_markers"

for (i in 1:nrow(pg_markers)){
  seurat_marker$celltype[seurat_marker$seurat_clusters %in% pg_markers$seurat_markers[i]] <- pg_markers$Consensus[i]
}
seurat_marker <- subset(seurat_marker, celltype != 'Unknown' & celltype != "Accessory Gland?" & celltype != "Muscle")

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

library(slingshot)
library(tradeSeq)


Xgenes <- filter(ortholog_table, chr == "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  c() %>% intersect(rownames(seurat_marker))

DefaultAssay(seurat_marker) <- "RNA"
Xexp <- PercentageFeatureSet(seurat_marker, features = Xgenes)
seurat_marker <- AddMetaData(object = seurat_marker, metadata = Xexp, col.name = 'xprop')

FeaturePlot(seurat_marker, features = 'xprop', split.by = 'treatment')


sce <- as.SingleCellExperiment(seurat_marker, assay = 'RNA')
sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC_and_Spermatogonia',
                 end.clus = 'Spermatids')
counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)
library(grDevices)
colorspt <- colorRampPalette(brewer.pal(11,'RdBu')[-6])(100)
plotcolpt <- colorspt[cut(sce$slingPseudotime_1, breaks=100)]

colorsct <- brewer.pal(8,'Set2')
names(colorsct) <- unique(sce$celltype)
colorsct2 <- brewer.pal(16,'Set3')
names(colorsct2) <- unique(sce$integrated_snn_res.0.2)
#use a brewer.pall pallete with 16 colours for colorsct2



plotcolct <- colorsct[sce$celltype]
plotcolct2 <- colorsct2[sce$integrated_snn_res.0.1]
dev.off()
layout(matrix(1:3, nrow = 1))
par(mar = c(4.5,4,1,1))
plot(reducedDims(sce)$UMAP, col = plotcolct, pch=16, asp = 1)
plot(reducedDims(sce)$UMAP, col = plotcolct2, pch=16, asp = 1)
plot(reducedDims(sce)$UMAP, col = plotcolpt, pch=16, asp = 1)
lines(SlingshotDataSet(sce)@curves$Lineage2, lwd=2, col='black')

seurat_marker <- AddMetaData(seurat_marker, pseudotime[,2], col.name = 'pseudotime')
plots <- list()
plots[[1]] <- seurat_marker@meta.data %>% 
  ggplot(aes(x = pseudotime, y = xprop, colour = treatment)) + #geom_point(alpha = 0.1) + 
  geom_smooth() + 
  theme_classic() + 
  scale_colour_brewer('Greens')

plots[[2]] <- seurat_marker@meta.data %>% 
  ggplot(aes(x = pseudotime, colour = treatment)) + geom_density() + 
  scale_colour_brewer('Greens') + 
  theme_classic()

plots[[3]] <- ggplot(seurat_marker@meta.data, aes(x = pseudotime, y = celltype)) + 
  geom_boxplot(outlier.alpha = 0)

ggarrange(plotlist = plots, nrow = 3, common.legend = TRUE, legend = "bottom")




figure <- ggarrange(plots[[1]]  + rremove("xlab"), plots[[2]] + rremove("xlab"), plots[[3]] + rremove("xlab"),# remove axis labels from plots
                    labels = NULL,
                    nrow = 3,
                    common.legend = TRUE, legend = "bottom",
                    align = "hv")


figure

ggplot(sce@metadata, aes(x = slingPseudotime_1, y = xprop, color = treatment)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'loess', se = FALSE, span = 0.1) +
  scale_color_manual(values = c('red', 'blue')) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = 'Pseudotime', y = 'X-linked gene expression') +
  facet_wrap(~celltype, ncol = 4)


############# TRADESEQ -------------
gois <- lapply(cell_type_DEG_output, function(x)(return(rbind(x$ST_bias, 
                                                          x$SR_bias)))) %>% 
  bind_rows()


sce <- fitGAM(sce)
# test for dynamic expression
icMat <- evaluateK(sce,
                   conditions = factor(sce$treatment),
                   nGenes = 100,
                   k = 4:5, parallel=F)

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
