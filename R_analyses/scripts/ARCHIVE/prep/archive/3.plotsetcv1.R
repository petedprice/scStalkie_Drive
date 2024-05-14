#!/usr/bin/env Rscript

###### SETTING UP INPUT COMMANDS ----
library("optparse")
option_list = list(
  make_option(c("-d", "--path_to_integrated_seurat_object"), type="character", default=".", 
              help="path to where you have the stored your seurat object", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default=".", 
              help="where you want to save your output plots and RData files", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="number of threads for parallelising", metavar="numeric"),
  make_option(c("-u", "--orthologs"), type="character", default=, 
              help="path to orthologs RData", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path_to_integrated_seurat_object)){
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
load(opt$path_to_integrated_seurat_object) # path to filtered seurat RData
load("outdata/RData/orthologs.RData")
markers <- filter(orthologs_testis, is.na(Cluster) == FALSE) %>% 
  filter(sub("^gene-", "", TDel_GID) %in% rownames(seurat_integrated))

#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)

#parallelise
plan("multicore", workers = opt$threads)
options(future.globals.maxSize = 8000 * 1024^5)

#PCAS/UMAPS/ETC?ETC
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:30,
                             reduction = "pca")

seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Visualization
#p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "treatment")
#p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE)
#ggsave("plots/UMAP1.pdf", p2 + p1, width = 15, height = 10)

DefaultAssay(seurat_integrated) <- "RNA"
#nk.markers <- FindConservedMarkers(seurat_integrated, ident.1 = 6, grouping.var = "treatment", verbose = FALSE)
#head(nk.markers)
########################### ORTHOLOGS ---------






plot_func <- function(cluster, mk_df = markers){
  print(cluster)
  mks <- filter(mk_df, Cluster == cluster)
  mks2 <- str_split(mks$TDel_GID, "gene-", simplify = TRUE)[,2]
  size = length(mks2) * 1.5
  f <- FeaturePlot(seurat_integrated, features = mks2, min.cutoff = "q10")
  ggsave(paste("plots/", cluster, "_feature_plot.pdf", sep = ""), f, height = size, width = 1.5* size)
}

lapply(unique(markers$Cluster), plot_func, mk_df = markers)


ggsave("plots/featureplot_markers.pdf", plots, height = 30, width = 30)


f <- FeaturePlot(seurat_integrated, features = marker_genes[51:80], min.cutoff = "q9")
ggsave("plots/featureplot_markers.pdf", f, height = 30, width = 30)
d <- DimPlot(seurat_integrated, reduction = "umap", split.by = "treatment", height = 25, width = 25)
ggsave("plots/del.pdf", d)
