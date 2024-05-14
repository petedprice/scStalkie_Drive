###LIBRARIES ----
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

#### FUNCTINOS ----
collapse <- function(x){
  gene <- t(str_split(x[15], ",", simplify = TRUE))
  gene <- gsub(" ", "", gene, fixed = TRUE)
  gene <- c(str_split(gene, ":", simplify = TRUE))
  gene <- gene[gene != ""]
  Cluster = rep(x[2], length(gene))
  return(data.frame())
}

plot_func<- function(cluster, mk_df = markers){
  print(cluster)
  mks <- filter(mk_df, Cluster == cluster)
  mks2 <- str_split(mks$TDel_GID, "gene-", simplify = TRUE)[,2]
  size = length(mks2) * 1.5
  plots <- lapply(mks2, FeaturePlot, object = seurat_integrated, min.cutoff = "q10")  
  pdf(paste(plotpath, cluster, "_feature_plot.pdf", sep = ""), width = 5, height = 5)
  for (p in plots){
    plot(p)
  }
  dev.off()
}
#### LOAD DATA ----
load("outdata/RData/integrated_seurat.RData")
markers_in <- readxl::read_excel("indata/markers/WITTPLOS2021/WITT_PLOS_GEN_2021_MARKERSv1.xlsx", col_names = FALSE)
markers_in <- readxl::read_excel("indata/markers/WITTPLOS2021/Compiled_markers13012022.xlsx", col_names = TRUE) %>% 
  pivot_longer(c("Marker1", "Marker2"), names_to = 'gene')
markers_in <- markers_in[,c("Cell_type", "value")]
markers_in <- markers_in[is.na(markers_in$value) == FALSE,]
load("outdata/RData/orthologs.RData")
##### SCRIPTS ----

## RUN seurat prep, including PCA, UMAP and clustering ----
#DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.4)

DimPlot(seurat_integrated, group.by = "ident")

### Read in markers ----
markers2df <- apply(markers_in, 1, function(x)(return(
  data.frame(
    Cluster = rep(x[1]), 
    gene = t(str_split(x[2], ",", simplify = TRUE)) %>% 
      gsub(pattern = " ", replacement =  "")) 
))) %>% bind_rows()

markers2df[markers2df == "P-cup"] <- 'p-cup'
markers <- merge(ortholog_table, markers2df, by.x = 'Dros_GID', by.y = 'gene', all.x = TRUE, all.y = TRUE)
markers = markers[is.na(markers$Cluster) == FALSE, ]

mks <- gsub("gene-", "", markers$TDel_GID)

for (m in mks){
  m2 <-gsub("gene-", "", m)
  if (m2 %in% rownames(seurat_integrated)){
    pdf(paste("plots/del/", m, "_del2.pdf", sep = ""))
    p <- FeaturePlot(object = seurat_integrated, features = m2, min.cutoff = "q10")
    plot(p)
    dev.off()
    stop()
  }
}

pdf("plots/del.pdf")
FeaturePlot(seurat_integrated, features = mks)
dev.off()

plotpath = "plots/"
lapply(unique(markers$Cluster), plot_func, mk_df = markers)
