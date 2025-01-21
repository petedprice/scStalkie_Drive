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

##### LOADING DATA ----
load(args[1])

print_barcode <- function(x){
  bc <- str_split(x, "_", simplify = T) %>% .[,2]
  write.table(bc, 
              paste0(x, "_cluster.txt"),
              row.names = F, col.names = F, quote = F)

}

lapply(unique(seurat_integrated$cells), print_barcode)
