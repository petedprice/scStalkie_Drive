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

# split object into a list by sample
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")


if (samples != 'all'){
  print("sample removing")
  keep_samples <- read.table(opt$samples)[,1]
  split_seurat <- split_seurat[keep_samples]
  print(paste("keeping samples ", keep_samples, sep = ""))
}

print(paste("memory available:", mem_used()))
if (doublet_finder == TRUE){
  ##### DoubletFinder ----- 
  cat("You chose to remove doublets so removing them for you!\n \
      To keep doublets use command -z FALSE when running this script")
  nExp <- lapply(split_seurat, function(x)return(round(ncol(x)*0.025))) 
  split_seurat <- lapply(split_seurat, SCTransform, vars.to.regress =
                                c("mitoRatio","nUMI","S.Score","G2M.Score"))

  print("SCT trnasformed")
  print(paste("memory available:", mem_used()))
  split_seurat <- lapply(split_seurat, RunPCA)
  split_seurat <- lapply(split_seurat, function(x)(return(RunUMAP(x, dims = 1:40,reduction = "pca"))))
  print("PCA and SCT done (1)")  
  split_seurat <- lapply(split_seurat, function(x)(return(FindNeighbors(x, dims = 1:40, verbose = FALSE))))
  split_seurat <- lapply(split_seurat, function(x)(return(FindClusters(x, verbose = FALSE))))
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
  
  #optimal_pks <- lapply(split_seurat, find_pk)
  optimal_pks <- list()
  for (i in 1:length(split_seurat)){
	print(i)
	optimal_pks[[i]] <-find_pk(split_seurat[[i]])
	print(i)
        print(paste("pregc memory available:", mem_used()))
	Sys.meminfo()
	gc()
	print(paste("postgc memory available:", mem_used()))
  	Sys.meminfo()
  }
  names(optimal_pks) <-names(split_seurat)
  print(optimal_pks)
  print("optimal pk found")
  split_seurat <- lapply(names(split_seurat), 
                         function(x)(doubletFinder_v3(split_seurat[[x]] , 
                                                      pN=0.25, pK=optimal_pks[[x]], nExp=nExp[[x]], PCs=1:40, sct=TRUE)))
  names(split_seurat) <-names(optimal_pks)  

  print("doubletfinder run")
  df_names <- lapply(split_seurat, function(x)(return(colnames(x@meta.data)[grepl("DF.classification", colnames(x@meta.data))])))
  print(df_names)
  
  ch_df_name <- function(so, df_name){
    so@meta.data <- so@meta.data %>% 
      dplyr::rename(doublet_finder = df_name)
    return(so)
  }
  
  split_seurat <- sapply(names(split_seurat), function(x)(return(ch_df_name(split_seurat[[x]], df_names[[x]]))))
  print("name change done")
  print(split_seurat)
  plot_doublets <- function(x){
    DimPlot(x, group.by = 'doublet_finder')
  }
  doublet_plots <- lapply(names(split_seurat), function(x)(return(DimPlot(split_seurat[[x]], group.by = 'doublet_finder') +
                                                                   ggtitle(x))))
  names(doublet_plots) <- names(optimal_pks)
  print("plots made")
  for (p in names(doublet_plots)){
	ggsave(filename = paste(plotpath, p, "_doublet_plots.pdf", sep = ""),plot = doublet_plots[[p]],  width = 10, height = 17)
  }
  print("plots saved")
  split_seurat <- lapply(split_seurat, function(x)(return(subset(x,doublet_finder == "Singlet"))))
  print("seur obj subsetted")
}

#SCT normalize the data (SCTransform also accounts for sequencing depth, 
#also regressing out mitochondrial percentage and cell cycle scoring)
print("SCTransform")
split_seurat <- lapply(split_seurat, SCTransform, vars.to.regress = 
                                c("mitoRatio","nUMI","S.Score","G2M.Score"))

gc()
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
print(paste("memory available:", mem_used()))
print("integrating")
#seurat_integrated <- IntegrateData(anchorset = anchors, 
#                                   normalization.method = "SCT", 
#                                   features.to.integrate = unique(unlist(lapply(split_seurat, rownames))))
print(paste("memory available:", mem_used()))
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

