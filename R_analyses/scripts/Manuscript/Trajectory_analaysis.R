library(Seurat)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
#remotes::install_version("matrixStats", version="1.1.0")
library(tradeSeq)
library(slingshot)

load("indata/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned.RData")
seurat_integrated$treatment <- "SR"
seurat_integrated$treatment[grep("st", seurat_integrated$sample)] <- "ST"
DefaultAssay(seurat_integrated) <- "integrated"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
celltype_table <- read.table("data/celltype_table.txt", header = T)
#using integrated_snn_res.0.4 from seurat_integrated, assign celltypes using celltype_table
seurat_integrated$celltype <- NA
for (i in 1:nrow(celltype_table)){
  seurat_integrated@meta.data$celltype[seurat_integrated@meta.data$integrated_snn_res.0.4 == celltype_table$clus_no[i]] <- celltype_table$celltype[i]
}


DimPlot(seurat_integrated, group.by = 'celltype')


#seurat_integrated <- subset(seurat_integrated, celltype != "Muscle" & celltype != "Epithelial")

DefaultAssay(seurat_integrated) <- "SCT"
Xgenes <- filter(ortholog_table, chr == "Chr_X") %>% 
  dplyr::select(REF_GENE_NAME) %>% unlist() %>% 
  gsub("_", "-", .) %>% 
  c() %>% intersect(rownames(seurat_integrated))

Xexp <- PercentageFeatureSet(seurat_integrated, features = Xgenes)
seurat_integrated <- AddMetaData(object = seurat_integrated, metadata = Xexp, col.name = 'xprop')

FeaturePlot(seurat_integrated, features = 'xprop', split.by = 'treatment')

sce <- as.SingleCellExperiment(seurat_integrated, assay = 'RNA')
sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC/Spermatogonia',
                 end.clus = 'Spermatids')
counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)
library(grDevices)
colorspt <- colorRampPalette(brewer.pal(11,'RdBu')[-6])(100)
plotcolpt <- colorspt[cut(sce$slingPseudotime_3, breaks=100)]

colorsct <- brewer.pal(8,'Set2')
names(colorsct) <- unique(sce$celltype)
plotcolct <- colorsct[sce$celltype]
dev.off()
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))
plot(reducedDims(sce)$UMAP, col = plotcolct, pch=16, asp = 1)
plot(reducedDims(sce)$UMAP, col = plotcolpt, pch=16, asp = 1)
lines(SlingshotDataSet(sce)@curves$Lineage3  , lwd=2, col='black')

seurat_integrated <- AddMetaData(seurat_integrated, pseudotime[,3], col.name = 'pseudotime')
plots <- list()
plots[[1]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = pseudotime, y = nFeature_SCT, colour = treatment)) + #geom_point(alpha = 0.1) + 
  geom_smooth() + 
  theme_classic() + 
  scale_colour_brewer('Greens')

plots[[2]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = pseudotime, y = xprop, colour = treatment)) + #geom_point(alpha = 0.1) + 
  geom_smooth() + 
  theme_classic() + 
  scale_colour_brewer('Greens')

plots[[3]] <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = pseudotime, colour = treatment)) + geom_density() + 
  scale_colour_brewer('Greens') + 
  theme_classic()

plots[[4]] <- ggplot(seurat_integrated@meta.data, aes(x = pseudotime, y = celltype)) + 
  geom_boxplot(outlier.alpha = 0)

ggarrange(plotlist = plots, nrow = 4, common.legend = TRUE, legend = "bottom")




figure <- ggarrange(plots[[1]]  + rremove("xlab"), plots[[2]] + rremove("xlab"), plots[[3]] + rremove("xlab"),# remove axis labels from plots
                    labels = NULL,
                    nrow = 3,
                    common.legend = TRUE, legend = "bottom",
                    align = "hv")

