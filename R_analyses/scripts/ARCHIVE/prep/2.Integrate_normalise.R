#!/usr/bin/env Rscript
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
library(DoubletFinder)
library("optparse")


###### SETTING UP INPUT COMMANDS ----
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
              help="path to dataframe containing cell cycle markers", metavar="character"),
  make_option(c("-z", "--doublet_finder"), type="character", default="TRUE", 
              help="whether to remove or keep doublets", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path_to_seurat_object)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



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


if (opt$doublet_finder == TRUE){
  ##### DoubletFinder ----- 
  cat("You chose to remove doublets so removing them for you!\n \
      To keep doublets use command -z FALSE when running this script")
  nExp <- future_lapply(split_seurat, function(x)return(round(ncol(x)*0.025))) 
  split_seurat <- future_lapply(split_seurat, SCTransform)
  split_seurat <- future_lapply(split_seurat, RunPCA)
  split_seurat <- future_lapply(split_seurat, function(x)(return(RunUMAP(x, dims = 1:40,reduction = "pca"))))
  
  split_seurat <- future_lapply(split_seurat, function(x)(return(FindNeighbors(x, dims = 1:40, verbose = FALSE))))
  split_seurat <- future_lapply(split_seurat, function(x)(return(FindClusters(x, verbose = FALSE))))
  
  find_pk <- function(so){
    sweep.list <- paramSweep_v3(so, PCs = 1:40, sct = T)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    return(optimal.pk)
  }
  
  optimal_pks <- future_lapply(split_seurat, find_pk)
  
  split_seurat <- sapply(names(split_seurat), 
                         function(x)(doubletFinder_v3(split_seurat[[x]] , 
                                                      pN=0.25, pK=optimal_pks[[x]], nExp=nExp[[x]], PCs=1:40, sct=TRUE)))
  
  df_names <- lapply(split_seurat, function(x)(return(colnames(x@meta.data)[grepl("DF.classification", colnames(x@meta.data))])))
  
  ch_df_name <- function(so, df_name){
    so@meta.data <- so@meta.data %>% 
      rename(doublet_finder = df_name)
    return(so)
  }
  
  split_seurat <- sapply(names(split_seurat), function(x)(return(ch_df_name(split_seurat[[x]], df_names[[x]]))))
  
  
  plot_doublets <- function(x){
    DimPlot(x, group.by = 'doublet_finder')
  }
  doublet_plots <- lapply(names(split_seurat), function(x)return(DimPlot(split_seurat[[x]], group.by = 'doublet_finder') +
                                                                   ggtitle(x)))
  
  dp2 <- ggarrange(plotlist = doublet_plots, nrow = 4, ncol = 2, common.legend = T)
  ggsave(filename = paste(plotpath, "doublet_plots.pdf", sep = ""),  width = 10, height = 17)
}

split_seurat <- future_lapply(split_seurat, function(x)(return(subset(x,doublet_finder == "Singlet"))))

#SCT normalize the data (SCTransform also accounts for sequencing depth, 
#also regressing out mitochondrial percentage and cell cycle scoring)
print("SCTransform")
split_seurat <- future_lapply(split_seurat, SCTransform, vars.to.regress = 
                                c("mitoRatio","nUMI","S.Score","G2M.Score"))


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

seurat_SCT_normalised <- FindNeighbors(seurat_SCT_normalised, dims = 1:30, verbose = FALSE)
seurat_SCT_normalised <- FindClusters(seurat_SCT_normalised, verbose = FALSE)

save(seurat_SCT_normalised, file = paste(outdatapath, "/seurat_SCT_normalised.RData", sep = ""))

fs_PCA1 <- DimPlot(seurat_SCT_normalised,
                   split.by = "sample")
ggsave(filename = paste(plotpath, "fs_sample_PCA.pdf", sep = ""))
fs_PCA2 <- DimPlot(seurat_SCT_normalised,
                   split.by = "treatment")
ggsave(filename = paste(plotpath, "fs_treatment_PCA.pdf", sep = ""))
ggsave(filename = paste(plotpath, "fs_treatment_PCA.pdf", sep = ""))
fs_PCA3 <- DimPlot(seurat_SCT_normalised,
                   split.by = "Phase", group.by = "Phase")
ggsave(filename = paste(plotpath, "fs_Phase_PCA.pdf", sep = ""))



#Integrate data 
print("integrating")
seurat_integrated <- IntegrateData(anchorset = anchors, 
                                   normalization.method = "SCT", 
                                   features.to.integrate = unique(unlist(lapply(split_seurat, rownames))))
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunICA(object = seurat_integrated)

seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
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

