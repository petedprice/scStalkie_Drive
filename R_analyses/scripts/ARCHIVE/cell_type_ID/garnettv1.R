###LIBRARIES ----
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(ggpubr)
library(devtools)
library(SeuratWrappers)
library(cicero)

# First install Bioconductor and Monocle 3
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install()

# Next install a few more dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment'))

devtools::install_github('cole-trapnell-lab/monocle3')
# Install a few Garnett dependencies:
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db'))
# Install Garnett
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")

library(garnett)

load("data/RData/seurat_markers.RData")
load("outdata/RData/ortholog_table.RData")

cds <- as.cell_data_set(seurat_integrated)
cds <- cluster_cells(cds, resolution=1e-3)
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
patchwork::wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
p1 <- plot_cells(cds,
           color_cells_by = "customclassif",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
p2 <- plot_cells(cds,
                 color_cells_by = "scina_labels",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)



pdf("plots/Cell_types/garnet_trajectory.pdf", width = 11, height = 4)
ggarrange(plotlist = list(p1,p2))
dev.off()

swap_names <- function(x, tab, srt){
  names <- unlist(lapply(x, function(g)(return(tab$TDel_GID[tab$Dros_GID == g])))) %>% 
    gsub(pattern = "gene-", replacement = "")
  return(intersect(names, rownames(srt)))
}
markers <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx", col_names = TRUE) %>% 
  dplyr::select("Gene", "Cluster")
markerslist <- lapply(unique(markers$Cluster), function(x)(return(markers$Gene[markers$Cluster == x])))
names(markerslist) <- unique(markers$Cluster)
names(markerslist) <- unique(markers$Cluster)
nmarkerslist <- sapply(markerslist, swap_names, tab = ortholog_table, srt = seurat_integrated)

garnett_table_func <- function(cell, ml){
  name <- paste(">", cell, sep = "")
  name <- gsub(",", " or ", name)
  genes <- paste(ml[[cell]], collapse = ", ")
  genes_exp <- paste("expressed: ", genes, sep = "")
  out_table <- matrix(data = c(name, genes_exp, " "))
  return(out_table)
}

garnett_table_matrix <- lapply(names(nmarkerslist), garnett_table_func, ml = nmarkerslist) %>% 
  unlist() %>% 
  as.data.frame()
write.table(garnett_table_matrix, file = "outdata/cell_ID/garnett_marker_table.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)


marker_file_path <- "outdata/cell_ID/garnett_marker_table.txt" 

marker_check <- check_markers(cds, marker_file_path,
                              db='none')

plot_markers(marker_check)

## training the classifier ----
cds_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = marker_file_path, 
                                        db = 'none', num_unknown= 50)


cds <- classify_cells(cds, cds_classifier,
                           db = 'none',
                           cluster_extend = F)
table(pData(cds)$cell_type)

