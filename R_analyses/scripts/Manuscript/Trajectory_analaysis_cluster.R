library(Seurat)
library(dplyr)
library(princurve)
library(TrajectoryUtils)
library(slingshot, lib = '.')
library(tradeSeq, lib = '.')
library(DelayedMatrixStats, lib = '.')

#install.packages("princurve", lib = '.')
#BiocManager::install("TrajectoryUtils", lib = '.')

args <- commandArgs(trailingOnly=TRUE)
load(args[1])
run_type = args[2]


Idents(seurat_final) <- seurat_final$celltype
seurat_final <- JoinLayers(seurat_final)
DefaultAssay(seurat_final) <- "RNA"
Idents(seurat_final) <- seurat_final$celltype
keep_clusters <- c("GSC/Spermatogonia", "Primary spermatocytes", "Spermatocytes", "Secondary spermatocytes", "Early spermatids", "Late spermatids", "Spermatids")
keep_clusters <- keep_clusters[keep_clusters %in% unique(seurat_final$celltype)]

# 
# rm_cells = c("sr5_GAGGGTACACCTAAAC-1", "sr5_GATAGAAAGCCTATCA-1",
# "st1_AAGTTCGGTATTGGCT-1", "sr5_AAGCATCTCGTGGCTG-1", "sr5_ACGTAGTAGTCAACAA-1", "sr5_CAAGGGATCTATTGTC-1", 
# "sr5_GACCCTTAGTGTCATC-1", "sr5_GACCGTGCATTGCTTT-1", "sr5_GATTCTTCACCAAATC-1", 
# "sr5_GCATCTCGTACACTCA-1", "sr5_GGACGTCGTAGGAGGG-1", "sr5_TAAGTCGAGACATAAC-1",
# "sr5_TCCAGAAAGGCAGGTT-1", "sr5_TCGCACTCACGCCACA-1", "sr5_TGTTGGACACAGACGA-1", "sr5_TTCGATTCAAAGGGTC-1")

sce <- seurat_final %>% 
  subset(., idents = keep_clusters) %>% 
  as.SingleCellExperiment(., assay = 'RNA')

#sce <- sce[,!colnames(sce) %in% rm_cells]

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
              cellWeights = cellWeights, nknots = 8, conditions = as.factor(sce$treatment))
  save(sce, sce_trad, file = "sce_GAMed_8k.RData")
} 




#################################################





