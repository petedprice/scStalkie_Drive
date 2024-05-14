#!/usr/bin/env Rscript

###LIBRARIES ----
library(Seurat)
library(tidyr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

##### Data and custom scripts ----
path_to_MD2022 <- args[1]
datapath = args[2]
output_path = args[3]
metadata_nf=read.csv(args[4], header = F)
mt_genes=read.table(args[5])
files <- list.files(datapath)
samples <- metadata_nf[,1]
files <- files[files %in% samples]
colnames(metadata_nf) <- c("sample", "species", "ref", "sex", "stage")
metadata_nf$files <- files

print(metadata_nf)

source(paste(path_to_MD2022, "/Rscripts/seurat/Usefull_functions.R", sep = ""))

#filtering thresholds 
filt = "filtered"
ftr = 200

#READING IN COUNT MATRICES and creating sample variables. 
c=0
for (file in metadata_nf$files){
  print(file)
  c=c+1
  seurat_data <- Read10X(data.dir = paste(datapath, "/", file, "/outs/filtered_feature_bc_matrix/", sep = ""))

  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = ftr, ### FILTERING CELLS WITH LESS THAN 200 GENES
                                   project = metadata_nf$sample[c], 
                                   min.cells = 3) 
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = metadata_nf$sample[c])
  seurat_obj@meta.data$sample <- metadata_nf$sample[c]
  if (c == 1){
	merged_seurat = seurat_obj
  } else {
	merged_seurat = merge(x = seurat_obj, y = merged_seurat)
  }
}

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA) # Calculating the number of features per UMI 




#CALCULATING PROPORTION OF READS MAPPING TO MITOCHONDRIAL GENOME (x100 as originally a percentage)

mt2=intersect(mt_genes$V1, rownames(merged_seurat))

merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, features = mt2)/ 100 
metadata <- merged_seurat@meta.data # Seperately creating metadata dataframe to save metrics etc without risking affecting seurat object
metadata$cells <- rownames(metadata) #add cell IDs to metadat
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
merged_seurat@meta.data <- metadata #Save the more complete metadat to the seurat object

#FILTERING DATA
filtered_seurat <- subset(x = merged_seurat, 
                           # (log10GenesPerUMI > 0.8) & # Can be dying cells or simple cell types such as blood cells
                            (mitoRatio < 0.05))
metadata_clean <- filtered_seurat@meta.data

outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
save(filtered_seurat, metadata_clean, file = paste(outdatapath, "/filtered_seurat.RData", sep = ""))

plotpath = paste(output_path, "/plots", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)
make_plots_function(metadata_clean, plotpath = plotpath)




