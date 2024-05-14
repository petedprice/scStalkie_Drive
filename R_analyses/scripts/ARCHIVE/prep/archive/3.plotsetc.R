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
  make_option(c("-u", "--orthologs"), type="character", default=NA, 
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
library(dplyr)

#Load data and parsing commands
output_path <- opt$output_path
load(opt$path_to_integrated_seurat_object) # path to filtered seurat RData
load(opt$orthologs)
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


plot_func <- function(cluster, mk_df = markers){
  print(cluster)
  mks <- filter(mk_df, Cluster == cluster)
  mks2 <- str_split(mks$TDel_GID, "gene-", simplify = TRUE)[,2]
  size = length(mks2) * 1.5
  plots <- lapply(mks2, FeaturePlot, object = seurat_integrated, min.cutoff = "q10")  
  pdf(paste(plotpath, cluster, "_feature_plot.pdf", sep = ""), width = 5, height = 5)
  for (p in plots){
    plot(p)
  }
  dev.off()
}
lapply(unique(markers$Cluster), plot_func, mk_df = markers)
plots <- list()
plots[[1]] <- DimPlot(seurat_integrated, group.by = "sample", pt.size = 1)
plots[[2]] <- DimPlot(seurat_integrated, split.by = "sample", ncol = 3)
plots[[3]] <- DimPlot(seurat_integrated, group.by = "ident", pt.size = 1)
save(plots,  file = paste(outdatapath, "/3.plots.RData", sep = ""))

pdf(paste(plotpath, "DimPlotst.pdf", sep = ""), width = 20, height = 20)
plot(plots[[1]])
plot(plots[[2]])
plot(plots[[3]])
dev.off()

save(seurat_integrated,  file = paste(outdatapath, "/3.plot_seurat.RData", sep = ""))

