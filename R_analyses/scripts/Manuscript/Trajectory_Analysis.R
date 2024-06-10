library(Seurat)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(slingshot)
library(tradeSeq)
library(grDevices)
library(TSCAN)
library(pheatmap)

load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_ss.RData")
load("data/trajectory/sce_GAMed.RData")
sce_trad <- sce
################################################

rowData(sce_trad)$assocRes <- associationTest(sce_trad, lineages = T, l2fc = log2(0.5))
assocRes <- rowData(sce_trad)$assocRes
assocRes <- assocRes[is.na(assocRes$waldStat) == F,]

assocRes %>% 
  ggplot(aes(x = pvalue_lineage1_conditionST, y = pvalue_lineage1_conditionSR)) + 
  geom_point()

bias_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionSR, "fdr") <= 0.05)
]

yhatSmooth <- predictSmooth(sce_trad, gene = c("gene-9128"), nPoints = 50, tidy = FALSE)
yhatSmooth <- yhatSmooth[is.na(yhatSmooth[,1]) == FALSE,]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = T, 
                       show_colnames = T)


condRes <- conditionTest(sce, l2fc = log2(2))
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

################# del
plotSmoothers(sce_trad, assays(sce)$counts,
              gene = 'PB.370',
              alpha = 1, border = TRUE)
################

# least significant gene
plotSmoothers(sce_trad, assays(sce_trad)$counts,
              gene = rownames(assays(sce_trad)$counts)[oo[nrow(sce_trad)]],
              alpha = 1, border = TRUE)


#################################################
Muscle <- c(17)
Spermatocytes <- c(10,11,13)
Spermatids <- c(3,4)
`GSC/Spermatogonia` <- c(0,2,7)
Pre_meiotic_cyst <- c(5,12,14,15)
Post_meiotic_cyst <- c(1,6,8,9,16)

### Plotting key markers against numbered cell clusters for assignment ---

seurat_integrated_ss@meta.data$celltype <- 'NA'
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Muscle] <- "Muscle"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatocytes] <- "Spermatocytes"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatids] <- "Spermatids"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% `GSC/Spermatogonia`] <- "GSC & Spermatogonia"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Pre_meiotic_cyst] <- "Pre-meiotic \ncyst"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Post_meiotic_cyst] <- "Post-meiotic \ncyst"
seurat_integrated_ss$treatment <- "SR"
seurat_integrated_ss$treatment[grep("st", seurat_integrated_ss$sample)] <- "ST"
Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype

DefaultAssay(seurat_integrated_ss) <- "integrated"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]


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
