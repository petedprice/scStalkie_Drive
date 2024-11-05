# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
load("data/RData/integrated_seurat.RData")
load("outdata/RData/ortholog_table.RData")


#### FUNCTIONS ----
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

DEG_func <- function(c, seurat_obj, ortholog_table){
  ss_obj <- subset(seurat_obj, customclassif == c)
  t1 <- unique(seurat_obj$treatment)[1]
  t2 <- unique(seurat_obj$treatment)[2]
  DEG_test <- FindMarkers(ss_obj, ident.1 = t1, 
                          ident.2 = t2, 
                          logfc.threshold = log(2))
  DEG_test$DEG_comp <- c
  DEG_test$gene <- rownames(DEG_test)
  rownames(DEG_test) <- NULL
  return(DEG_test)
}


markers <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx", col_names = TRUE) %>% 
  dplyr::select("Gene", "Cluster")

## PREPPING DATA ----
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")


seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.4)


markerslist <- lapply(unique(markers$Cluster), function(x)(return(markers$Gene[markers$Cluster == x])))
names(markerslist) <- unique(markers$Cluster)
nmarkerslist <- sapply(markerslist, swap_names, tab = ortholog_table, srt = seurat_integrated)
gs_list <- list()
gs_list$gs_positive <- nmarkerslist


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


seurat_integrated@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_integrated@meta.data$customclassif[seurat_integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#### SCINA ----
keep_or_not <- lapply(nmarkerslist, check_subset, ml = nmarkerslist) %>% 
  unlist()
SCINA_markerslist <- nmarkerslist[keep_or_not]


scina.data <- as.data.frame(seurat_integrated@assays$integrated[,]) 

results = SCINA(scina.data, SCINA_markerslist, 
                max_iter = 1, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
seurat_integrated$scina_labels <- results$cell_labels


## PLOTS ----
d1 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE, split.by = 'treatment')        
d2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE, 
             group.by = 'customclassif', split.by = 'treatment')  
d3 <- DimPlot(seurat_integrated, reduction = 'umap', group.by = "scina_labels", label = T, split.by = 'treatment')
d <- ggarrange(plotlist = list(d1,d2,d3), nrow = 3)
pdf("plots/Cell_types/SCINA_sc_type.pdf", width = 22, height = 24)
d
dev.off()

test_seurat <- subset(seurat_integrated, nUMI > 20000)

Idents(seurat_integrated) <- seurat_integrated$treatment
### ID differentially expressed genes! ----
DEG_results <- future.apply::future_lapply(unique(test_seurat$customclassif), DEG_func, 
                             seurat_obj = test_seurat,
                             ortholog_table = ortholog_table) %>% 
  bind_rows()

DEG_results$ref_gene <- paste("gene-", DEG_results$gene, sep = "")
DEG_results <- DEG_results %>% left_join(ortholog_table[,c(3,5,6)], by = c('ref_gene' = 'TDel_GID'))
DEG_results %>% filter(p_val_adj < 0.05) %>% 
  ggplot(aes(fill = Dros_CHR, x = DEG_comp)) + geom_bar()


DEG_results %>% filter(p_val_adj < 0.05) %>% 
  group_by(DEG_comp) %>% 
  summarise(n = n())
 
