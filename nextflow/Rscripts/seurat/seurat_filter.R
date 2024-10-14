#!/usr/bin/env Rscript
###LIBRARIES ----
library(Seurat)
library(tidyr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

##### Data and custom scripts ----
path_to_pd <- args[1]
datapath = args[2]
output_path = args[3]
metadata_nf=read.csv(args[4], header = F)
mt_genes=read.table(args[5])
nfs=as.numeric(args[6])
mtr=as.numeric(args[7])
gu=as.numeric(args[8])
mt_genes$V1 <- gsub("_", "-", mt_genes$V1)
print(metadata_nf)
samples <- metadata_nf[,1]
print(samples)
colnames(metadata_nf) <- c("sample", "species", "ref", "sex", "type", "mt_contig")

print(metadata_nf)

source(paste(path_to_pd, "/Rscripts/seurat/Usefull_functions.R", sep = ""))

#filtering thresholds 
filt = "filtered"
ftr = nfs

#READING IN COUNT MATRICES and creating sample variables. 
seurat_data <- Read10X(data.dir = datapath)

seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = ftr, ### FILTERING CELLS WITH LESS THAN 200 GENES
                                   project = metadata_nf$sample, 
                                   min.cells = 3)
 
seurat_obj <- RenameCells(seurat_obj, add.cell.id = metadata_nf$sample)
seurat_obj@meta.data$sample <- metadata_nf$sample
merged_seurat = seurat_obj

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA) # Calculating the number of features per UMI 




#CALCULATING PROPORTION OF READS MAPPING TO MITOCHONDRIAL GENOME (x100 as originally a percentage)

mt2=intersect(mt_genes$V1, rownames(merged_seurat))

merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, features = mt2)/ 100 
metadata <- merged_seurat@meta.data # Seperately creating metadata dataframe to save metrics etc without risking affecting seurat object
metadata$mitofilter=mtr
metadata$featurefilter=nfs
metadata$complexityfilter=gu
metadata$cells <- rownames(metadata) #add cell IDs to metadat
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
merged_seurat@meta.data <- metadata #Save the more complete metadat to the seurat object

#FILTERING DATA
filtered_seurat <- subset(x = merged_seurat, 
                            (log10GenesPerUMI > gu) & # Can be dying cells or simple cell types such as blood cells
                            (mitoRatio < mtr)) #& 
#				(nGene < 7000))
metadata_clean <- filtered_seurat@meta.data

outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
save(filtered_seurat, metadata_clean, file = paste(outdatapath, "/filtered_seurat.RData", sep = ""))

plotpath = paste(output_path, "/plots", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)
make_plots_function(metadata_clean, plotpath = plotpath)




