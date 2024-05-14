## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(ggpubr)
#library(scater)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library("biomaRt")
library(clusterProfiler)
library(ggrepel)
library(lme4)
library(GenomicFeatures)

load("data/RData/marker_seurat.RData")
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
seurat_marker$celltype = seurat_marker$sctype_labels
seurat_marker$treatment = "ST" #which samples start with SR treatment is SR
seurat_marker$treatment[grepl("sr", seurat_marker$sample)] = "SR"
#DimPlot(seurat_marker,  label = TRUE, group.by = 'celltype')


#SC_TYPE FUNCTIONS
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


##### LOADING DATA ----


colums <- colnames(seurat_marker@meta.data)
clusters <- seurat_marker@meta.data[,colums[startsWith(colums, "inte")]] %>%
  as.data.frame() %>%
  apply(., 2, function(x)(return(max(as.numeric(x)))))

resolution <- which.min(abs(clusters - 24)) %>% names()


seurat_marker$seurat_clusters <- seurat_marker[[resolution]]
markers <- read.csv("outdata/markers_elife.csv")

clusters <- unique(markers[,'celltype'])
clusters <- clusters[is.na(clusters) == F]

markerslist <- lapply(clusters, function(x)(return(markers$marker[markers[,'celltype'] == x])))
markerslist <- lapply(markerslist, function(x)(return(x[which(is.na(x) == F)])))
names(markerslist) <- clusters

gs_list <- list()
gs_list$gs_positive <- markerslist

#### SC_TYPE MARKERS ----
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seurat_marker[["integrated"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat_marker@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_marker@meta.data[seurat_marker@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_marker@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[])

seurat_marker@meta.data$sctype_labels = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  seurat_marker@meta.data$sctype_labels[seurat_marker@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seurat_marker, group.by = 'sctype_labels', label = T)
seurat_marker$celltype <- seurat_marker$sctype_labels
