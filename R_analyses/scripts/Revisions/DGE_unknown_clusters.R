################# load libraries and functions and data prep ###################
library(dplyr)
library(Seurat)
library(openxlsx)
library(ggpubr)
library(scuttle)
library(edgeR)
library(statmod)
library(scran)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(ggrepel)
library(lme4)
library(GenomicFeatures)
library(ggrepel)
library(lme4)
rm(list = ls())

load("data/RData/seurat_unfiltered.RData")
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated$treatment <- "SR"

seurat_integrated$keep <- seurat_integrated$celltype
seurat_integrated$celltype <- seurat_integrated$celltype2

seurat_integrated$treatment[grep("st", seurat_integrated$sample)] <- "ST"
ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table <- ortholog_table %>% group_by(REF_GENE_NAME, OF_DMEL, 
                                              FBgnOF, OMA_REFGENE, chr, OMA_CG, 
                                              OMA_DMEL, consensus_gene) %>% 
  summarise(start = min(start), end = max(end))
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_integrated)) %>% unique()
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  gsub("_", "-", .) %>% 
  intersect(rownames(seurat_integrated))%>% unique()

seurat_integrated <- JoinLayers(seurat_integrated)
sce <- seurat_integrated %>% as.SingleCellExperiment(assay = "RNA")

################################################################################


######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]


## NEW PBDGE FUNCTION
PBDGE <- function(label, sce, func = 'custom', seurat_obj = seurat_integrated, consistent_genes = T){
  if (consistent_genes == TRUE){
    ST_cells <- filter(seurat_integrated@meta.data, celltype == label & treatment == "ST") %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    SR_cells <- filter(seurat_integrated@meta.data, celltype == label & treatment == "SR") %>% 
      dplyr::select(cells) %>% c() %>% unlist()
    if (length(ST_cells) < 2 | length(SR_cells) < 2){
      run = F
    } else {
      run = T
      ST_keep <- (rowSums(seurat_integrated@assays$RNA$counts[,ST_cells] > 1) > (0.05 * length(ST_cells))) #| 
      #rowSums(seurat_integrated@assays$RNA$counts[,ST_cells] > 1) > 10)
      
      SR_keep <- (rowSums(seurat_integrated@assays$RNA$counts[,SR_cells] > 1) > (0.05 * length(SR_cells)))# | 
      #                  rowSums(seurat_integrated@assays$RNA$counts[,SR_cells] > 1) > 10)
      
      keep <- ST_keep | SR_keep
    }
  }else {
    keep <- rep(T, nrow(sce))
  }
  if (run == F){
    return("Not enough cells in one of the groups")
  }
  dge_data <- sce[keep,sce$celltype == label]
  y <- DGEList(counts(dge_data), samples=colData(dge_data), group = dge_data$treatment)
  y <- calcNormFactors(y)
  design <- model.matrix(~treatment, y$samples)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  cpm <- cpm(y, log = F)
  if (func == "custom"){
    keep1 <- cpm[,y$samples$treatment== "ST"] >= 1
    keep1 <- rowSums(keep1) > (ncol(keep1)/2)
    keep2 <- cpm[,y$samples$treatment== "SR"] >= 1
    keep2 <- rowSums(keep2) > (ncol(keep2)/2)
    keep <- keep1 | keep2
  } else if (func == "edger"){  
    keep <- filterByExpr(y, group = y$sample$treatment, 
                         large.n = 3, min.prop = 0.6, min.count = 5)
    #keep <- filterByExpr(y, group = y$sample$treatment)
  }
  y <- y[keep,keep.lib.sizes=FALSE]
  
  STSR <- exactTest(y, pair = c("ST", "SR")) 
  fit <- glmQLFit(y, design, robust = T)
  res <- glmQLFTest(fit, coef=ncol(design))
  out <- list()
  out$res <- res
  out$logcpm <- cpm(y, log = T, prior.count = .0001)
  out$cpm <- cpm(y, log = F)
  out$toptags <- topTags(res, n = nrow(res))
  out$exacttest <- topTags(STSR, n = nrow(STSR))
  out$func <- func
  out$y <- y
  return(out)
}


for (ct in unique(agg_cell_types$celltype)){
  print(ct)
  PBDGE(ct, sce = agg_cell_types, func = 'custom', seurat_obj = seurat_integrated, 
        consistent_genes = T)
}



de.results_act_edger <- sapply(unique(agg_cell_types$celltype), PBDGE, sce = agg_cell_types, func = 'custom',
                               USE.NAMES = T, simplify = F)

table(seurat_integrated$celltype, seurat_integrated$sample)
cols <- c("treatment", "celltype", "sample")
get_cpm_values_func <- function(x){
  if (!is.character(x)){
    md <- x$res
    ST <- rownames(filter(md$samples, treatment == "ST"))
    SR <- rownames(filter(md$samples, treatment == "SR"))
    cpm <- x$logcpm %>% 
      as.data.frame() %>%
      rownames_to_column("genes") %>% 
      pivot_longer(cols = colnames(.)[2:ncol(.)], names_to = "sample", values_to = "logcpm") %>% 
      mutate(Chr = ifelse(genes %in% Xgenes, "X", 
                          ifelse(genes %in% Agenes, "A", "Mito"))) %>% 
      mutate(treatment = ifelse(sample %in% ST, "ST", "SR")) %>% 
      mutate(celltype = md$samples$celltype[1]) %>% 
      rename(edgeRsample = sample) %>% 
      mutate(sample = md$samples[edgeRsample,]$sample)
    return(cpm)
  }
}  


tidy_cpm_edger <- lapply(de.results_act_edger, get_cpm_values_func) %>% bind_rows()
tidy_cpm_edger$func <- "custom"

################ Medians approach for dosage compensation ------
medianAs <- tidy_cpm_edger %>% 
  group_by(celltype, treatment, genes, Chr) %>%
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  group_by(celltype, treatment) %>%
  summarise(medianA = median(logcpm[Chr == "A" & logcpm >= 1]), 
            n = length(logcpm[Chr == "A" & logcpm >= 1]))


XA <- tidy_cpm_edger %>% 
  merge(medianAs) %>% 
  group_by(celltype, genes, Chr, treatment, medianA) %>% 
  summarise(logcpm = log2(mean(2^logcpm))) %>% 
  #group_by(sample, celltype, treatment, func, Chr) %>%
  summarise(XA = logcpm - medianA, 
            logcpm = logcpm )



## No of DGE ##


########################  ENRICHMENT FIGURE ##########################
get_DGE_data <- function(ct,x){
  if (!is.character(x[[ct]])){
    out = x[[ct]][['toptags']]$table
    out$celltype <-  ct
    out$genes <- rownames(out)
    rownames(out) <- NULL
    return(out)
  }
}

dif_exp_data <- lapply(names(de.results_act_edger), get_DGE_data, x = de.results_act_edger) %>% 
  bind_rows() %>% 
  mutate(Significant = ifelse(.$logFC > 1 & .$FDR < 0.05, "ST-biased",
                              ifelse(.$logFC < -1 & .$FDR < 0.05, "SR-biased", "Unbiased"))) %>% 
  merge(ortholog_table, by.x = 'genes', by.y = 'REF_GENE_NAME')

unknown_DGE <- dif_exp_data %>% 
  filter(celltype %in% c("Unknown 1", "Unknown 2", "Unknown 4")) %>% 
  filter(Significant != "Unbiased") %>% 
  dplyr::select(c(genes, logFC, FDR, celltype, Significant, consensus_gene))

write.table(unknown_DGE, "data/Unknown_celltypes_DGE.csv", 
            sep = ",", row.names = F, quote = F)


dif_exp_data_table <-  dif_exp_data %>% 
  mutate(chr = ifelse(chr %in% c("Chr_1", "Chr_2"), "Autosome", "X")) %>% 
  group_by(celltype, Significant, chr) %>% 
  mutate(celltype = gsub("\n", "", celltype)) %>% 
  summarise(count = n()) %>% 
  reshape2::dcast(celltype + chr ~ Significant, value.var = "count") %>% 
  rename(`Cell type` = celltype, Chromosome = chr)


dif_exp_data_table[is.na(dif_exp_data_table)] <- 0


a <- seurat_integrated %>% DimPlot(., split.by = 'treatment', group.by = 'celltype', label = T)
b <- seurat_integrated@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA, fill = treatment)) + 
  geom_boxplot()
ggarrange(a,b, nrow = 2)

c <-table(seurat_integrated$celltype, seurat_integrated$sample) %>% 
  ggtexttable()
d <- dif_exp_data_table %>% 
  ggtexttable()

          
##############################################################





####################
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]

#"Unknown 5"
#"Late spermatids",
agg_cell_types_new <- agg_cell_types[,agg_cell_types$celltype %in% c("Unknown 4",  "Early spermatids") & agg_cell_types$treatment == "SR"]
agg_cell_types_new$celltype[agg_cell_types_new$celltype %in% c("Unknown 5", "Unknown 4")] <- "Unknown"
agg_cell_types_new$celltype[agg_cell_types_new$celltype %in% c("Late spermatids", "Early spermatids")] <- "Spermatids"
## NEW PBDGE FUNCTION

dge_data <- agg_cell_types_new
y <- DGEList(counts(dge_data), samples=colData(dge_data), group = dge_data$celltype)
y <- calcNormFactors(y)
design <- model.matrix(~celltype, y$samples)
y <- normLibSizes(y)
y <- estimateDisp(y, design)
cpm <- cpm(y, log = F)

keep1 <- cpm[,y$samples$celltype== "Spermatids"] >= 1
keep1 <- rowSums(keep1) > (ncol(keep1)/2)
keep2 <- cpm[,y$samples$celltype== "Unknown"] >= 1
keep2 <- rowSums(keep2) > (ncol(keep2)/2)
keep <- keep1 | keep2

y <- y[keep,keep.lib.sizes=FALSE]

STSR <- exactTest(y, pair = c("Spermatids", "Unknown")) 
fit <- glmQLFit(y, design, robust = T)
res <- glmQLFTest(fit, coef=ncol(design))
out <- list()
out$res <- res
out$logcpm <- cpm(y, log = T, prior.count = .0001)
out$cpm <- cpm(y, log = F)
out$toptags <- topTags(res, n = nrow(res))
out$exacttest <- topTags(STSR, n = nrow(STSR))
out$y <- y


out$res$table %>% 
  mutate(gene = rownames(.)) %>% 
  merge(ortholog_table, by.x = 'gene', by.y = "REF_GENE_NAME") %>% 
  View()

