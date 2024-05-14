#!/usr/bin/env Rscript

###### SETTING UP INPUT COMMANDS ----
library("optparse")
option_list = list(
  make_option(c("-d", "--path_to_seurat_object"), type="character", default=".", 
              help="path to where you have the stored your seurat object", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default=".", 
              help="where you want to save your output plots and RData files", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="number of threads for parallelising", metavar="numeric"),
  make_option(c("-s", "--samples"), type="character", default="all", 
              help="path to dataframe containing samples (see format on github)", metavar="character"),
  make_option(c("-c", "--cellcycle"), type="character", default="/home/bop20pp/software/MeioticDrive2022/R_analyses/data/cell_cycle_markers_complete.csv", 
              help="path to dataframe containing cell cycle markers", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path_to_seurat_object)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)

library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(future)
library(future.apply)

#Load data and parsing commands
output_path <- opt$output_path
load(opt$path_to_seurat_object) # path to filtered seurat RData
cellcycle <- read.table(opt$cellcycle)
cellcycle$TDel_GID <- gsub("gene-", "", cellcycle$TDel_GID)

#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)

#parallelise
plan("multicore", workers = opt$threads)
options(future.globals.maxSize = 8000 * 1024^5)

#Cell cycle scoring 
filtered_seurat <- CellCycleScoring(filtered_seurat, 
                                    g2m.features = cellcycle$TDel_GID[cellcycle$phase == "G2/M"], 
                                    s.features = cellcycle$TDel_GID[cellcycle$phase == "S"])

# split object into a list by sample
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

if (opt$samples != 'all'){
  print("sample removing")
  keep_samples <- read.table(opt$samples)[,1]
  split_seurat <- split_seurat[keep_samples]
  print(paste("keeping samples ", keep_samples, sep = ""))
}

#SCT normalize the data (SCTransform also accounts for sequencing depth, 
#also regressing out mitochondrial percentage and cell cycle scoring)
print("SCTransform")
split_seurat <- future_lapply(split_seurat, SCTransform, vars.to.regress = 
                                c("mitoRatio","nUMI","S.Score","G2M.Score"))
ss_names <- c('normal', 'st', 'sr')
split_seurat_st <- split_seurat[startsWith(names(split_seurat), "st")]
split_seurat_sr <- split_seurat[startsWith(names(split_seurat), "sr")]
splits <- list(split_seurat, split_seurat_st, split_seurat_sr)
names(splits) <- ss_names
names(ss_names) <- ss_names

#prep data for integration 
# Identify variable features for integrating
print("SelectIntegrationFeatures")
features <- future_lapply(splits, SelectIntegrationFeatures, nfeatures = 3000, 
                   USE.NAMES = TRUE)

#Prepossessing step necessary if SCT transformed
print("PrepSCTIntegration")
splits <- future_lapply(ss_names, function(x)(return(PrepSCTIntegration(object.list = splits[[x]],  
                                                                 anchor.features = features[[x]]))))

#Find anchors that link datasets
print("FindIntegrationAnchors")
anchors <- future_lapply(ss_names, function(x)(return(FindIntegrationAnchors(object.list = splits[[x]], 
                                                                      anchor.features = features[[x]], 
                                                                      normalization.method = "SCT"))))

### INTEGRATING 
print("integrating")
integrateds <- future_lapply(ss_names, function(x)(return(IntegrateData(anchorset = anchors[[x]], 
                                                                 normalization.method = "SCT", 
                                                                 features.to.integrate = 
                                                                   unique(unlist(lapply(splits[[x]], rownames)))))))
seurat_SCT_normaliseds <- future_lapply(ss_names, function(x)(return(Reduce(merge, splits[[x]]))))

umap_tsne_pca <- function(x, f = NULL){
  x <- RunPCA(object = x, features = f)
  x <- RunTSNE(x, dims = 1:40, reduction = "pca")
  x <- RunUMAP(x,  dims = 1:40, reduction = "pca")
  x <- FindNeighbors(object = x, dims = 1:40)
  x <- FindClusters(object = x, resolution = 0.4)
}

integrateds <- future_lapply(integrateds, umap_tsne_pca)
seurat_SCT_normaliseds <- future_lapply(ss_names, function(x)(return(umap_tsne_pca(seurat_SCT_normaliseds[[x]], f = features[[x]]
                                                                            ))))

save(anchors, features, seurat_SCT_normaliseds, file = paste(outdatapath, "/features_anchors_split.RData", sep = ""))
save(seurat_SCT_normaliseds, file = paste(outdatapath, "/SCT_normaliseds_split.RData", sep = ""))
save(integrateds,file = paste(outdatapath, "/integrateds_split.RData", sep = ""))

