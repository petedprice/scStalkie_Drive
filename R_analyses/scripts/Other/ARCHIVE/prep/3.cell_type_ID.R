#!/usr/bin/env Rscript
##### LIBRARIES -------
library("optparse")
library(SCINA)
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

#SC_TYPE FUNCTIONS
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

###### SETTING UP INPUT COMMANDS ----
option_list = list(
  make_option(c("-d", "--path_to_seurat_object"), type="character", default=".", 
              help="path to where you have the stored your seurat object", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default=".", 
              help="where you want to save your output plots and RData files", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="number of threads for parallelising", metavar="numeric"),
  make_option(c("-l", "--ortholog_table"), type="character", default="/home/bop20pp/software/MeioticDrive2022/R_analyses/data/ortholog_table.txt", 
              help="path to dataframe containing ortholog information", metavar="character"),
  make_option(c("-s", "--marker_source"), type="character", default="comp_clusters", 
              help="which markers to use for classifying cell types. Column name from ortholog_table.", metavar="character")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
output_path <- opt$output_path

if (is.null(opt$path_to_seurat_object)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

#make folders
outdatapath = paste(output_path, "/outdata", sep = "")
dir.create(outdatapath, showWarnings = F, recursive = T)
plotpath = paste(output_path, "/plots/", sep = "")
dir.create(plotpath, showWarnings = F, recursive = T)

#parallelise
plan("multicore", workers = opt$threads)
options(future.globals.maxSize = 8000 * 1024^5)

##### FUNCTIONS ----
swap_names <- function(x, tab, srt){
  names <- unlist(lapply(x, function(g)(return(tab$TDel_GID[tab$Dros_GID == g])))) %>% 
    gsub(pattern = "gene-", replacement = "")
  return(intersect(names, rownames(srt)))
}

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
load(opt$path_to_seurat_object)
seurat_integrated$seurat_clusters <- seurat_integrated$integrated_snn_res.0.1
ortholog_table <- read.table(opt$ortholog_table)
marker_source <- opt$marker_source
clusters <- unique(ortholog_table[,marker_source])
clusters <- clusters[is.na(clusters) == F]

markerslist <- lapply(clusters, function(x)(return(ortholog_table$TDel_GID[ortholog_table[,marker_source] == x])))
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

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

seurat_integrated@meta.data$sctype_labels = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_integrated@meta.data$sctype_labels[seurat_integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


#### SCINA ----
keep_or_not <- lapply(markerslist, check_subset, ml = markerslist) %>% 
  unlist()
SCINA_markerslist <- markerslist[keep_or_not]
SCINA_markerslist <- markerslist

scina.data <- as.data.frame(seurat_integrated@assays$integrated[,]) 

results = SCINA(scina.data, SCINA_markerslist, 
                max_iter = 1, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
seurat_integrated$scina_labels <- results$cell_labels


seurat_marker <- seurat_integrated
seurat_marker@meta.data$marker_source <- marker_source
d1 <- DimPlot(seurat_marker, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'sctype_labels')  +
  ggtitle("SC_type clusters")+  theme(plot.title = element_text(hjust = 0.5))

d2 <- DimPlot(seurat_marker, reduction = 'umap', group.by = "scina_labels", label = T) +
  ggtitle("SCINA clusters")+  theme(plot.title = element_text(hjust = 0.5))
d <- ggarrange(plotlist = list(d1,d2), nrow = 1)
ggsave(filename = paste(plotpath, "CELL_TYPES.pdf", sep = ""),  width = 17, height = 8.5)

d1 <- DimPlot(seurat_marker, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'sctype_labels', split.by = 'treatment')  +
  ggtitle("SC_type clusters")+  theme(plot.title = element_text(hjust = 0.5))
d2 <- DimPlot(seurat_marker, reduction = 'umap', group.by = "scina_labels", label = T, 
              split.by = 'treatment') +
  ggtitle("SCINA clusters")+  theme(plot.title = element_text(hjust = 0.5))
d <- ggarrange(plotlist = list(d1,d2), nrow = 2)
ggsave(filename = paste(plotpath, "CELL_TYPES_treatment_split.pdf", sep = ""),  width = 17, height = 17)

save(seurat_marker, file = paste(outdatapath, "/", marker_source, "_marker_seurat.RData", sep = ""))

