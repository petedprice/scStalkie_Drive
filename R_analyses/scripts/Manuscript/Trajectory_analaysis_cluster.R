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
eval_k = args[2]

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

DefaultAssay(seurat_integrated_ss) <- "RNA"
Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype


sce <- seurat_integrated_ss %>% 
  subset(., idents = c("GSC & Spermatogonia", "Spermatocytes", "Spermatids")) %>% 
  as.SingleCellExperiment(., assay = 'RNA')

sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC & Spermatogonia',
                 end.clus = 'Spermatids')
print("check 3")

counts <- counts(sce)[which(rowSums(counts(sce)) != 0),]
print("check 4")
pseudotime <- slingPseudotime(sce, na = FALSE)
print("check 5")
cellWeights <- slingCurveWeights(sce)
print("check 6")

################################# TRADESEQ -------------------------------------
#evaluate k 
if (eval_k == "Yes"){
  print(eval_k)
  icMat <- evaluateK(sce,
                     conditions = factor(sce$treatment),
                     nGenes = 200,
                     k = 3:10, parallel=F)
  
  save(icMat, "icMat.RData")
} else {
  icMat <- readRDS("icMat.RData")
  #FItting GAM
  sce <- fitGAM(counts = counts(sce), pseudotime = pseudotime, 
              cellWeights = cellWeights, nknots = args[3], conditions = as.factor(sce$treatment))
  save(sce, "sce_GAMed.RData")
}




#################################################





