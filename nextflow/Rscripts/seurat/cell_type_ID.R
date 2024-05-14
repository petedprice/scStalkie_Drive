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

#install.packages(c("SCINA"), repos = 'http://cran.us.r-project.org', lib = '.')
#library(SCINA,lib.loc = '.')
#SC_TYPE FUNCTIONS
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

###### SETTING UP INPUT COMMANDS ----
args = commandArgs(trailingOnly=TRUE)
outdatapath = paste(args[2], "/outdata/", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(args[2], "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)

##### FUNCTIONS ----

check_subset <- function(cell, ml){
  matches <- lapply(ml, function(x)(return(prod(cell %in% x)))) %>% 
    unlist()
  if (sum(matches) > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}



##### LOADING DATA ----
load(args[1])

colums <- colnames(seurat_integrated@meta.data)
clusters <- seurat_integrated@meta.data[,colums[startsWith(colums, "inte")]] %>% 
  as.data.frame() %>% 
  apply(., 2, function(x)(return(max(as.numeric(x)))))

resolution <- which.min(abs(clusters - 24)) %>% names()


seurat_integrated$seurat_clusters <- seurat_integrated[[resolution]]
markers <- read.csv(args[3])
markers$marker <- gsub("_", "-", markers$marker)

clusters <- unique(markers[,'celltype'])
clusters <- clusters[is.na(clusters) == F]

markerslist <- lapply(clusters, function(x)(return(markers$marker[markers[,'celltype'] == x])))
markerslist <- lapply(markerslist, function(x)(return(x[which(is.na(x) == F)])))
names(markerslist) <- clusters

gs_list <- list()
gs_list$gs_positive <- markerslist

#### SC_TYPE MARKERS ----
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seurat_integrated[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat_integrated@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_integrated@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
seurat_marker <- seurat_integrated
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

seurat_integrated@meta.data$sctype_labels = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_integrated@meta.data$sctype_labels[seurat_integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


#### SCINA ----
#keep_or_not <- lapply(markerslist, check_subset, ml = markerslist) %>% 
#  unlist()
#SCINA_markerslist <- markerslist[keep_or_not]
#SCINA_markerslist <- markerslist
#
#scina.data <- as.data.frame(seurat_integrated@assays$integrated[,]) 
#
#results = SCINA(scina.data, SCINA_markerslist, 
#                max_iter = 1, convergence_n = 10, 
#                convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
#                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
#seurat_integrated$scina_labels <- results$cell_labels


seurat_marker <- seurat_integrated
d1 <- DimPlot(seurat_marker, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'sctype_labels')  +
  ggtitle("SC_type clusters")+  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste(plotpath, "SC_TYPE_CELL_TYPES.pdf", sep = ""),  d1, width = 17, height = 8.5)

#d2 <- DimPlot(seurat_marker, reduction = 'umap', group.by = "scina_labels", label = T) +
#  ggtitle("SCINA clusters")+  theme(plot.title = element_text(hjust = 0.5))
#ggsave(filename = paste(plotpath, "SCINA_CELL_TYPES.pdf", sep = ""),  d2, width = 17, height = 8.5)


save(seurat_marker, file = paste(outdatapath, "/marker_seurat.RData", sep = ""))

