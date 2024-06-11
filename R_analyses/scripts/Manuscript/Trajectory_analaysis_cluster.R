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
rm_spermatid_cells <- c("sr5_AACAGGGGTAATGCTC-1", "sr5_AATGAAGCAAATGCTC-1", "sr5_AATGCCACAATGACCT-1", "sr5_AGCCACGCATCGTGGC-1",
                        "sr5_ATCACAGTCCTTCAGC-1", "sr5_ATCATTCTCGTGCACG-1", "sr5_CCCTAACCAACGCATT-1", "sr5_CTCATGCCACTGGATT-1",
                        "sr5_GCGGAAACAGAACATA-1", "sr5_GGACGTCTCTCGACGG-1", "sr5_GTGCTTCCACACGCCA-1", "sr5_TATCGCCGTCTCCTGT-1",
                        "sr5_TGCATCCAGCGACTTT-1", "sr5_TGCTTGCTCCGTATGA-1", "sr5_TGTCAGAGTATCAGGG-1", "st1_GGGTGAATCCGTTGGG-1",
                        "st2_AGGCCACGTCGAATTC-1", "st3_AGGCCACTCTGCTTTA-1", "st3_ATGCGATCAAGCGCTC-1", "sr3_TTGCGTCAGCGGATCA-1")

seurat_integrated_ss <- subset(seurat_integrated_ss, cells = rm_spermatid_cells, invert = TRUE)


Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype

DefaultAssay(seurat_integrated_ss) <- "RNA"
Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype


sce <- seurat_integrated_ss %>% 
  subset(., idents = c("GSC & Spermatogonia", "Spermatocytes", "Spermatids")) %>% 
  as.SingleCellExperiment(., assay = 'RNA')

sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = "UMAP", start.clus = 'GSC & Spermatogonia',
                 end.clus = 'Spermatids')

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





