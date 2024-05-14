#!/usr/bin/env Rscript
###LIBRARIES ----
library(Seurat)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggplot2)
library(future)
library(parallel)

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
output_path <- args[1]
#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)


seurat_objs <- list.files(pattern = ".RData")
split_seurat <-list()
for (sobj in seurat_objs){
	load(sobj)
	print("SCT transform")
	#doublet_seurat <- SCTransform(doublet_seurat, vars.to.regress =
        #                        c("nUMI","S.Score","G2M.Score"))
	doublet_seurat <- SCTransform(doublet_seurat, vars.to.regress =
                                c("mitoRatio","nUMI","S.Score","G2M.Score"))
	split_seurat[[doublet_seurat$sample[1]]] <-doublet_seurat
	
}



#### #prep data for integration ---- 
# Identify variable features for integrating
print("SelectIntegrationFeatures")
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

#Prepossessing step necessary if SCT transformed
print("PrepSCTIntegration")
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = features)

#Find anchors that link datasets
print("FindIntegrationAnchors")
anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                  anchor.features = features, 
                                  normalization.method = "SCT")

#Initial plot making for comparisons to after integration 

##### first Perform PCA -----
print("first unintegrated UMAPS")
seurat_SCT_normalised <- Reduce(merge, split_seurat)
seurat_SCT_normalised <- RunPCA(object = seurat_SCT_normalised, features = features)
seurat_SCT_normalised <- RunUMAP(seurat_SCT_normalised, 
                                 dims = 1:30,
                                 reduction = "pca")
gc()
seurat_SCT_normalised <- FindNeighbors(seurat_SCT_normalised, dims = 1:30, verbose = FALSE)
seurat_SCT_normalised <- FindClusters(seurat_SCT_normalised, verbose = FALSE)

save(seurat_SCT_normalised, file = paste(outdatapath, "/seurat_SCT_normalised.RData", sep = ""))

fs_PCA1 <- DimPlot(seurat_SCT_normalised,
                   split.by = "sample")
ggsave(filename = paste(plotpath, "fs_sample_PCA.pdf", sep = ""))
fs_PCA2 <- DimPlot(seurat_SCT_normalised,
                   split.by = "Phase", group.by = "Phase")
ggsave(filename = paste(plotpath, "fs_Phase_PCA.pdf", sep = ""))


gc()
#Integrate data 
print("integrating")
#seurat_integrated <- IntegrateData(anchorset = anchors, 
#                                   normalization.method = "SCT", 
#                                   features.to.integrate = unique(unlist(lapply(split_seurat, rownames))))
seurat_integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT")

print("running PCA")
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

print("running UMAP")
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")


seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                                 0.75, 1, 1.25, 1.5, 1.75, 2, 
                                                 2.5, 3))
#Save data
print("saving data")
save(seurat_integrated, file = paste(outdatapath, "/integrated_seurat.RData", sep = ""))
save(split_seurat, anchors, features, file = paste(outdatapath, "/integrated_seurat_extras.RData", sep = ""))

