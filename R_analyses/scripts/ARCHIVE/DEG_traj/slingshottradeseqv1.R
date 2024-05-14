if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tradeSeq")
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
library(tidyverse)
library(devtools)
install_github('kstreet13/bioc2020trajectories')
library(pheatmap)

load("data/RData/testis.Cluster_marker_seurat.RData")
load("data/DEG.RData")
ortholog_table <- read.table("data/ortholog_table.txt")

# run cell-type IDing
#siss <- subset(seurat_integrated, customclassif != 'Unknown' &
#                 integrated_snn_res.0.4 %in% c('0',1,2,3,4,5,6,7,8,9,10,11))
siss <- subset(seurat_marker, customclassif != 'Unknown')
siss <- seurat_marker
siss$new_clusters <- siss$customclassif
siss$new_clusters[siss$customclassif %in% c("Mature spermatids", "Early spermatids", 
                                            "Early spermatocytes")] <- "Spermatids"
siss$new_clusters[siss$customclassif == "GSC, Early spermatogonia"] <- "GSC/Spermatagonia"

DimPlot(siss, group.by = 'new_clusters', reduction = 'umap', split.by = 'treatment')

sce_trad <- as.SingleCellExperiment(siss, assay = "RNA")

shuffle <- sample(ncol(sce_trad))
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce_trad)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:2)[factor(colData(sce_trad)$treatment)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:2, bty = "n", 
       legend = levels(factor(colData(sce_trad)$treatment)))
cols = 1:length(unique(colData(sce_trad)$new_clusters))
plot(reducedDims(sce_trad)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(cols[factor(colData(sce_trad)$new_clusters)][shuffle], alpha = .5))
legend("topright", pch = 16, col = cols, bty = "n", legend = levels(factor(colData(sce_trad)$new_clusters)))

scores <- bioc2020trajectories::imbalance_score(
  rd = reducedDims(sce_trad)$UMAP, 
  cl = colData(sce_trad)$treatment,
  k = 20, smooth = 40)


grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce_trad)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)


sce_trad <- slingshot(sce_trad, reducedDim = 'UMAP', clusterLabels = sce_trad$celltype,
                 start.clus = "Cyst", approx_points = 20)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_trad$slingPseudotime_1, breaks=100)]

plotcol <- sce_trad$celltype

#pdf("plots/SMBE2023/TRAJ_ATLAS.pdf")
plot(reducedDims(sce_trad)$UMAP)
lines(SlingshotDataSet(sce_trad), color = 'red')
#dev.off()

sum_data <- data.frame(cluster = sce_trad$new_clusters, pt = sce_trad$slingPseudotime_1, 
           Type = sce_trad$treatment, 
           `Cell Type` = sce_trad$new_clusters)

sum_data$pseudotime <- sum_data$pt

sum_data$Type[sum_data$Type == 'sr'] <- "Drive"
sum_data$Type[sum_data$Type == 'st'] <- "Standard"

sum_data %>% 
  ggplot(aes(x = pseudotime, fill = Cell.Type)) + geom_density(alpha = 0.2) + 
  facet_grid(.~pseudotime)

sum_data %>% 
  ggplot(aes(x = Pseudotime, fill = Type)) + geom_density(alpha = 0.2)

pseudotime <- slingPseudotime(sce_trad, na = FALSE)
cellWeights <- slingCurveWeights(sce_trad)
gois <- lapply(cell_type_DEG_output, function(x)(return(c(x$ST_bias_genes$geneID, 
                                                          x$SR_bias_genes$geneID)))) %>% 
  unlist() %>% unique()
markers <- ortholog_table$TDel_GID[is.na(ortholog_table$comp_clusters) == FALSE]
null <- sample(rownames(sce)[(rownames(sce) %in% gois) == F], length(gois))
gois_null <- c(gois, null, markers) %>% unique()
length(gois_null)
#gois <- cell_type_DEG_output$`Mature spermatids`$ST_bias_genes
#gois_null <- c('LOC119669221', 'LOC119683206', 'LOC119681754', 
#  'LOC119677649', 'LOC119677648')


### ASSESSING NUMBER OF KNOTS TO USE FOR MODEL 
icMat <- evaluateK(counts = as.matrix(assays(sce_trad)$counts),
                   pseudotime = pseudotime,
                   cellWeights = cellWeights,
                   conditions = factor(colData(sce_trad)$treatment),
                   nGenes = 100,
                   k = 4:5, parallel=F)

counts <- counts(sce_trad)[which(rowSums(counts(sce_trad)) != 0),]
gois_null <- intersect(rownames(counts), gois_null)
#### FITTING GAM MODEL 
sce_trad_ss <- fitGAM(counts = counts[gois_null,], 
              pseudotime = pseudotime, 
              cellWeights = cellWeights,
              conditions = factor(colData(sce_trad)$treatment),
              nknots = 5, parallel=F)



## Testing for differentially expressed genes as a function of pseudotime
rowData(sce_trad_ss)$assocRes <- associationTest(sce_trad_ss, lineages = T, l2fc = log2(0.5))
assocRes <- rowData(sce_trad_ss)$assocRes
assocRes <- assocRes[is.na(assocRes$waldStat) == F,]

assocRes %>% 
  ggplot(aes(x = pvalue_lineage1_conditionST, y = pvalue_lineage1_conditionSR)) + 
  geom_point()

sr_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionSR, "fdr") <= 0.05)
]
st_genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionST, "fdr") <= 0.05)
]

#plot unique and overlap of genes that are DE across pseudotime between treatments. 
UpSetR::upset(fromList(list(st = st_genes, sr = sr_genes)))



yhatSmooth <- predictSmooth(sce_trad_ss, gene = st_genes, nPoints = 50, tidy = FALSE)
yhatSmooth <- yhatSmooth[is.na(yhatSmooth[,1]) == FALSE,]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = T, 
                       show_colnames = T)
 

### TESTING FOR DIFFERENTIAL EXPRESSION BETWEEN TREATMENTS ----
condRes <- conditionTest(sce_trad_ss, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
plotSmoothers(sce_trad_ss, assays(sce_trad_ss)$counts,
              gene = rownames(assays(sce_trad_ss)$counts)[oo[1]],
              alpha = 1, border = TRUE)

# least significant gene
plotSmoothers(sce_trad_ss, assays(sce_trad_ss)$counts,
              gene = rownames(assays(sce_trad_ss)$counts)[oo[nrow(sce_trad_ss)]],
              alpha = 1, border = TRUE)


##### HEATMAPS OF GENES DE BETWEEN CONDITIONS -----
### based on mean smoother
yhatSmooth <- predictSmooth(sce_trad_ss, gene = conditionGenes, nPoints = 50, tidy = FALSE)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_ST <- pheatmap(yhatSmoothScaled[, 51:100],
                           cluster_cols = FALSE,
                           show_rownames = F, show_colnames = T, main = "Standard", legend = FALSE,
                           silent = TRUE
)

matchingHeatmap_SR <- pheatmap(yhatSmoothScaled[heatSmooth_ST$tree_row$order, 1:50],
                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = T, show_colnames = T, main = "Drive",
                                 legend = FALSE, silent = TRUE
)

grid.arrange(heatSmooth_ST[[4]], matchingHeatmap_SR[[4]], ncol = 2)




