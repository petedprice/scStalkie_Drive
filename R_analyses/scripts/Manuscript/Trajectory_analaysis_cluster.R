library(Seurat)
library(dplyr)
library(princurve, lib = '.')
library(TrajectoryUtils, lib = '.')

library(slingshot, lib = '.')
library(tradeSeq, lib = '.')
#install.packages("princurve", lib = '.')
#BiocManager::install("TrajectoryUtils", lib = '.')

args <- commandArgs(trailingOnly=TRUE)
load(args[1])
run_type = args[2]


Idents(seurat_final) <- seurat_final$celltype

DefaultAssay(seurat_final) <- "RNA"
seurat_final@meta.data$celltype[seurat_final$integrated_snn_res.0.4 == 5] <- "Late spermatids"
Idents(seurat_final) <- seurat_final$celltype
keep_clusters <- c("GSC/Spermatogonia", "Primary Spermatocytes", "Spermatocytes", "Secondary Spermatocytes", "Spermatids", "Late spermatids")
keep_clusters <- keep_clusters[keep_clusters %in% unique(seurat_final$celltype)]

sce <- seurat_final %>% 
  subset(., idents = keep_clusters) %>% 
  as.SingleCellExperiment(., assay = 'RNA')

sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC/Spermatogonia',
                 end.clus = 'Late spermatids')

counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)

################################# TRADESEQ -------------------------------------
#evaluate k 
if (run_type == "eval_k"){
  print('evaluating K')
  pdf("k_evaluation.pdf", height = 7, width = 28)
  icMat <- evaluateK(counts, pseudotime = pseudotime, cellWeights = cellWeights,
                     conditions = factor(sce$treatment),
                     nGenes = 200,
                     k = 3:20)
  dev.off()
  
  save(icMat, file = "icMat.RData")
} else if (run_type == "fitGAM") {
  icMat <- load(file = "icMat.RData")
  #FItting GAM
  sce_trad <- fitGAM(counts = counts(sce), pseudotime = pseudotime, 
              cellWeights = cellWeights, nknots = 10, conditions = as.factor(sce$treatment))
  save(sce, sce_trad, file = "sce_GAMed.RData")
} 




#################################################





