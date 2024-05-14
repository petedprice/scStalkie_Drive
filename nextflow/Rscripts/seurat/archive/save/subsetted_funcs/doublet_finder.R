#!/usr/bin/env Rscript
###LIBRARIES ----
library(Seurat)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggplot2)
install.packages(c("fields", "spam", "pryr", "memuse"), repos = 'http://cran.us.r-project.org', lib = '.')
library(spam, lib.loc = '.')
library(fields, lib.loc = '.')
library(pryr, lib.loc = '.')
library(memuse, lib.loc = '.')
library(future)
library(parallel)

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
load(args[1]) # path to filtered seurat RData
output_path <- args[2]
threads = as.numeric(args[3])-2
samples="all"
cellcycle <- read.table(args[4])
doublet_finder = args[5]
doublet_finder = T
memory=340

print(paste("memory available:", mem_used()))
Sys.meminfo()
detectCores()

options(future.globals.maxSize = (memory*1000) * 1024^2)
plan("multiprocess", workers = threads)

if (doublet_finder == TRUE){
   remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade = F, lib = '.')
   library(DoubletFinder,lib.loc  = '.')
}

sessionInfo() 

#Load data and parsing commands

#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)


#Cell cycle scoring 
filtered_seurat <- CellCycleScoring(filtered_seurat, 
                                    g2m.features = cellcycle$gene[cellcycle$phase == "G2/M"], 
                                    s.features = cellcycle$gene[cellcycle$phase == "S"])

print(paste("memory available:", mem_used()))
if (doublet_finder == TRUE){
  ##### DoubletFinder ----- 
  cat("You chose to remove doublets so removing them for you!\n \
      To keep doublets use command -z FALSE when running this script")
  nExp <- round(ncol(filtered_seurat)*0.025) 
  filtered_seurat <-SCTransform(filtered_seurat, vars.to.regress =
                                c("mitoRatio","nUMI","S.Score","G2M.Score"))
  print("SCT trnasformed")
  filtered_seurat <-RunPCA(filtered_seurat)
  filtered_seurat <-RunUMAP(filtered_seurat, dims = 1:40,reduction = "pca")
  filtered_seurat <-FindNeighbors(filtered_seurat, dims = 1:40,verbose = F)
  filtered_seurat <- FindClusters(filtered_seurat, verbose = F)
  print("clusters and neighbours found (1)")
  find_pk <- function(so){
    sweep.list <- paramSweep_v3(so, PCs = 1:40, sct = T, num.cores = threads)
    print("sl")
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    print("ss")
    bcmvn <- find.pK(sweep.stats)
    print("bcmvn")
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    print("bcmvnmax")
    optimal.pk <- bcmvn.max$pK
    print("oppk")
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    return(optimal.pk)
  }
  
  optimal_pks <-find_pk(filtered_seurat)
  print(optimal_pks)
  print("optimal pk found")
  filtered_seurat <-doubletFinder_v3(filtered_seurat, pN=0.25, pK=optimal_pks, nExp=nExp, PCs=1:40, sct=TRUE)

  print("doubletfinder run")

  df_name <- colnames(filtered_seurat@meta.data)[grepl("DF.classification", colnames(filtered_seurat@meta.data))]

  ch_df_name <- function(so, df_name){
    so@meta.data <- so@meta.data %>% 
      dplyr::rename(doublet_finder = df_name)
    return(so)
  }
  filtered_seurat <- ch_df_name(filtered_seurat, df_name)
  print("name change done")
  doublet_plots <-DimPlot(filtered_seurat, group.by = "doublet_finder") + ggtitle(filtered_seurat$sample[1])

  ggsave(filename = paste(plotpath, filtered_seurat$sample[1], "_doublet_plots.pdf", sep = ""),plot = doublet_plots,  width = 10, height = 17)
 
  print("plots saved")
  doublet_seurat <- subset(filtered_seurat, doublet_finder == "Singlet")
  print("seur obj subsetted")
}

save(doublet_seurat, file = paste(outdatapath, "/doublet_seurat.RData", sep = ""))

