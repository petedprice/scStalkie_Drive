# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(SCINA)
library(ggpubr)
library(scater)
library(scuttle)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
load("data/RData/integrated_seurat.RData")
load("outdata/RData/ortholog_table.RData")
ortholog_table <- read.table("data/ortholog_table.txt")

####FUNCTIONS ----
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

DEG_func <- function(c, seurat_obj, ortholog_table, s = 'customclassif'){
  if (s == "customclassif"){
    ss_obj <- subset(seurat_obj, customclassif == c)
  } else if (s == "scina_labels") {
    ss_obj <- subset(seurat_obj, scina_labels == c)
  }
  if (length(unique(ss_obj$treatment)) < 2){
    return(NULL)
  } else {
    t1 <- unique(seurat_obj$treatment)[1]
    t2 <- unique(seurat_obj$treatment)[2]
    DEG_test <- FindMarkers(ss_obj, ident.1 = t1, 
                            ident.2 = t2, 
                            logfc.threshold = log(2), 
                            min.pct = 0.5)
    DEG_test$DEG_comp <- c
    DEG_test$gene <- rownames(DEG_test)
    rownames(DEG_test) <- NULL
    return(DEG_test)
  }
}


markers <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx", col_names = TRUE) %>% 
  dplyr::select("Gene", "Cluster")


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
d1 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE, split.by = 'treatment') + 
  ggtitle("PCA clusters") +  theme(plot.title = element_text(hjust = 0.5))

d2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE, 
             group.by = 'customclassif')  +
  ggtitle("SC_type clusters")+  theme(plot.title = element_text(hjust = 0.5))

d3 <- DimPlot(seurat_integrated, reduction = 'umap', group.by = "scina_labels", label = T) +
  ggtitle("SCINA clusters")+  theme(plot.title = element_text(hjust = 0.5))
d <- ggarrange(plotlist = list(d2,d3), nrow = 1)

pdf("plots/Cell_types/SCINA_sc_type.pdf", width = 9, height =7.8)
d2
dev.off()

pdf('plots/Cell_types/UMAP_pooled.pdf', width = 11, height = 8)
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("PCA clusters") +  theme(plot.title = element_text(hjust = 0.5))
dev.off()



### ID differentially expressed genes! ----
Idents(seurat_integrated) <- seurat_integrated$treatment
DEG_results1 <- future.apply::future_lapply(unique(seurat_integrated$scina_labels), DEG_func, 
                             seurat_obj = seurat_integrated,
                             ortholog_table = ortholog_table, 
                             s = 'scina_labels') %>% 
  bind_rows()
DEG_results2 <- future.apply::future_lapply(unique(seurat_integrated$customclassif), DEG_func, 
                                            seurat_obj = seurat_integrated,
                                            ortholog_table = ortholog_table, 
                                            s = 'customclassif') %>% 
  bind_rows()

comp_DEG <- left_join(DEG_results1, DEG_results2, by = c('gene' = 'gene', "DEG_comp" = "DEG_comp"), 
                      suffix = c("_SCINA", "_SC"))

comp_DEG$ref_gene <- paste("gene-", comp_DEG$gene, sep = "")

comp_DEG <- comp_DEG %>% left_join(ortholog_table[,c(3:6)], by = c('ref_gene' = 'TDel_GID'))
write.table(comp_DEG, "outdata/cellwise_DEG.txt")


pdf("plots/Cell_types/DGE.pdf", width = 16, height = 5)
comp_DEG %>% filter(DEG_comp != 'unknown') %>% 
  ggplot(aes(x = avg_log2FC_SCINA,
             y  = -log10(p_val_adj_SCINA),
             colour = (abs(avg_log2FC_SCINA) > 2 &  p_val_adj_SCINA < 0.05))) + 
  geom_point() + 
  facet_grid(~ DEG_comp) +
  geom_text(data=subset(comp_DEG, -log10(p_val_adj_SCINA) > 25 &abs(avg_log2FC_SCINA) > 2),
            aes(avg_log2FC_SCINA,-log10(p_val_adj_SCINA),label=Dros_GID)) + 
  ylab("-log10 Pvalue") + 
  xlab("log2 fold-change") + 
  theme(legend.position = "none")      
dev.off()

save(seurat_integrated, file = "data/RData/seurat_markers.RData")
#### Identifying sperm expression Y and not expressing the Y 
Y_expression_score <- function(x){
  
}

X_linked_genes <- filter(TDel_txdf, TDel_CHR == 'NC_051848.1')$TDel_GID %>% 
  unique() %>% 
  intersect(rownames(seurat_integrated))

DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated$Xexp <- (PercentageFeatureSet(object = seurat_integrated, features = X_linked_genes)/100)
seurat_integrated$Yexp <- (PercentageFeatureSet(assay = "RNA", object = seurat_integrated, features = Y_linked_genes)/100)
DefaultAssay(seurat_integrated) <- "integrated"


pdf("plots/Cell_types/X_linked_expression.pdf", width = 24, height = 4)
FeaturePlot(seurat_integrated, features = 'Xexp', split.by = 'customclassif')
dev.off()

pdf("plots/Cell_types/X_linked_expression.pdf", width = 8, height = 8)
FeaturePlot(seurat_integrated, features = 'Xexp')
dev.off()


seurat_integrated$Ysperm <- FALSE
seurat_integrated$Ysperm[seurat_integrated$Xexp < 0.0001] <- TRUE

pdf("plots/Cell_types/Y_sperm.pdf", width = 8, height = 8)
DimPlot(seurat_integrated, group.by = 'Ysperm', split.by = 'treatment')
dev.off()

DimPlot(seurat_integrated, group.by = 'nUMI')


FeaturePlot(seurat_integrated, features = rownames(seurat_integrated)[1:100])


seurat_integrated$expression <- NA
seurat_integrated$expression[seurat_integrated$nCount_SCT < 2500] <- "LOW"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 2500 & 
                               seurat_integrated$nCount_SCT < 5000] <- "MEDLOW"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 5000 & 
                               seurat_integrated$nCount_SCT < 7500] <- "MED"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 7500 & 
                               seurat_integrated$nCount_SCT < 10000] <- "HIGH"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 10000 & 
                               seurat_integrated$nCount_SCT < 12500] <- "VHIGH"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 12500] <- "TOP"

seurat_integrated$expression <- NA
seurat_integrated$expression[seurat_integrated$nCount_SCT < 5000] <- "LOW"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 5000 & 
                               seurat_integrated$nCount_SCT < 12500] <- "MED"
seurat_integrated$expression[seurat_integrated$nCount_SCT > 12500] <- "TOP"

seurat_integrated$complexity = seurat_integrated$nGene/seurat_integrated$nUMI
DimPlot(seurat_integrated, group.by = 'ident')

seurat_integrated@meta.data <- seurat_integrated@meta.data %>% 
  mutate(cat = cut(complexity, 
                   c(-Inf, quantile(complexity, c(0.01, .1, .25, .5, .75, .99)), Inf), 
                   labels = c("Super Low", "Low", "lowish", "Medium", "High", "super high", "wow")))




##### BIOCONDUCTOR DGESEQ2 -----
library(edgeR)
library(statmod)

sce <- as.SingleCellExperiment(seurat_integrated, assay = "RNA")
summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("treatment", "sample")])
summed
label = "Cyst"
#current <- summed[,label==summed$
                    
                    
                    
              
y <- DGEList(counts(summed), samples=colData(summed))
discarded <- summed$ncells < 10
y <- y[,!discarded]
summary(discarded)

keep <- filterByExpr(y, group=summed$treatment)
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

par(mfrow=c(2,4))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}
dev.off()
plotMDS(cpm(y, log=TRUE), 
        col= as.numeric(factor(y$samples$treatment)))


design <- model.matrix(~factor(treatment), y$samples)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

