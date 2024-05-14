#!/usr/bin/env Rscript
##### LIBRARIES -------
library(Seurat)
library(tidyr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(dplyr)
library(ggplot2)
#library(ggpubr)

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
outdatapath = paste(args[2], "/outdata/", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
clus = args[3]

##### LOADING DATA ----
load(args[1])


samp <- unique(seurat_integrated$sample)

print_barcodes_func <- function(samp, clus){
  barcodes <- filter(seurat_integrated@meta.data, 
                     sample == samp)
  barcodes2 <- barcodes[,c("cells", clus)]
  for (i in unique(barcodes2[,clus])){
    bc_out <- barcodes2 %>% 
      filter(!!sym(clus) == i) %>% 
      rownames(cells) %>% 
     str_split("_", simplify = T) %>% 
      .[,2]

    write.table(bc_out, 
                paste0(samp, "_", i, "cluster.txt"), 
                row.names = F, col.names = F, quote = F)
  }
}

lapply(samp, print_barcodes_func, clus)

