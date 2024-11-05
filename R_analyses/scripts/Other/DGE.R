## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(SCINA)
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


load("indata/RData/marker_seurat.RData")
ortholog_table <- read.table("data/ortholog_table.txt")

new_names <- ortholog_table$Dros_GID %>% 
  setNames(gsub("gene-", "", ortholog_table$TDel_GID))
new_names[is.na(new_names) == TRUE] <- ortholog_table$product[is.na(new_names)==T]
new_names[new_names == "missing"] <- ortholog_table$TDel_GID[new_names == "missing"]

ortholog_table$TDel_CHR_Name <- "unknown"
ortholog_table$TDel_CHR_Name[ortholog_table$TDel_CHR == "NC_051846.1"] <- 1
ortholog_table$TDel_CHR_Name[ortholog_table$TDel_CHR == "NC_051847.1"] <- 2
ortholog_table$TDel_CHR_Name[ortholog_table$TDel_CHR == "NC_051848.1"] <- "X"


##### BIOCONDUCTOR EDGER2 -----
## PLOT FUNCTION ---- 
make_pca_function <- function(x, ortholog_table){
  md <- metadata(x)
  mds_data <- md$y %>% 
    cpm(log = TRUE) %>% 
    plotMDS()
  
  PCA <- data.frame(x = mds_data$x, 
                    y = mds_data$y, 
                    names = md$y$samples$seq_folder,
                    treatment = md$fit$samples$treatment) %>% 
    ggplot(aes(x = x, y = y, fill = names, label = names, colour = treatment))  +
    geom_text(hjust="inward", vjust="inward") +
    labs(x = "PC2", y = "PC1", title = md$y$samples$new_clusters[1]) 
  pdf(paste("plots/DEG/PCA_", md$y$samples$new_clusters[1], ".pdf", sep = ""))
  print(PCA) 
  dev.off()
  x <- x[is.na(x$logFC) == FALSE,] %>%
    as.data.frame()
  
  
  x$geneID <- rownames(x)
  x$gene <- do.call(recode, c(list(rownames(x)), new_names))
  x <- merge(x, ortholog_table[,c('TDel_GID', 'TDel_CHR_Name')], by.x = 'geneID', by.y = 'TDel_GID')   
  top_p <- x$FDR[order(x$FDR)][5]
  top_p2 <- x$FDR[order(x$FDR)][5]
  
  Volcano <- x %>% 
    ggplot(aes(x = logFC, y = -log10(FDR), label = gene,
               colour = (abs(logFC) > 1 &  FDR < 0.05))) + 
    geom_point() + geom_text_repel(data=subset(x, FDR < top_p & abs(logFC) > 1),
                                   aes(logFC,-log10(FDR),label=gene), 
                                   hjust="inward", vjust="inward") +
    labs(y = "-log10 Pvalue (FDR adjusted)", x = "log2 fold-change", 
         title = paste(md$y$samples$new_clusters[1], "(positive is ST bias)"))+ 
    theme(legend.position = "none")  + 
    geom_hline(yintercept = -log10(0.05), color = 'grey', linetype = 'dashed') + 
    geom_vline(xintercept = c(-1,1), color = 'grey', linetype = 'dashed')
  pdf(paste("plots/DEG/DEG_", md$y$samples$new_clusters[1], ".pdf", sep = ""))
  print(Volcano) 
  dev.off()
  ST_bias_genes <- filter(x, logFC > 1, FDR < 0.05)
  SR_bias_genes <- filter(x, logFC < -1, FDR < 0.05)
  output <- list(PCA,Volcano, ST_bias_genes, SR_bias_genes, x)
  names(output) <- c("PCA", "Volcano", "ST_bias", "SR_bias", "all")
  return(output)
  
}

sce <- subset(seurat_marker, customclassif != 'Unknown') %>% 
  as.SingleCellExperiment(assay = "RNA")
sce$new_clusters <- sce$customclassif
sce$new_clusters[sce$customclassif %in% c("Mature spermatids", "Early spermatids", 
                                          "Early spermatocytes")] <- "Germ"
sce$new_clusters[sce$customclassif == "GSC, Early spermatogonia"] <- "Germ"

#### BULK RNASEQ ACROSS WHOLE TISSUE ----
whole_tissue <- aggregateAcrossCells(sce, id=colData(sce)[,c("sample")])
de.results_wt <- pseudoBulkDGE(whole_tissue, 
                               label='bulk',
                               condition = whole_tissue$treatment,
                               design=~treatment,
                               coef = "treatmentst",
)
metadata(de.results_wt$bulk)$y$samples$sctype_labels <- "bulk"

bulk_DEG_output <- make_pca_function(de.results_wt$bulk, ortholog_table)
DEG <- bulk_DEG_output$all %>% 
  mutate(DGE)

### PSEUDOBULK INDIVIDUAL CELL TYPES ----
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("new_clusters", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]

y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
keep <- filterByExpr(y, group=agg_cell_types$treatment)
y <- y[keep,]
mds_data <-cpm(y, log = TRUE) %>% 
  plotMDS()
plot <- data.frame(x = mds_data$x, 
                   y = mds_data$y, 
                   names = y$samples$seq_folder,
                   treatment = y$samples$treatment, 
                   cell_type = y$samples$new_clusters) %>% 
  ggplot(aes(x = x, y = y, label = cell_type, colour = treatment))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "ind cell types")

pdf("plots/DEG/PCA_pseudobulk_ind_cell_types.pdf")
plot 
dev.off()

de.results_act <- pseudoBulkDGE(agg_cell_types, 
                                label=agg_cell_types$new_clusters,
                                condition = agg_cell_types$treatment,
                                design=~treatment,
                                coef = "treatmentst"
)



cell_type_DEG_output <- lapply(de.results_act, make_pca_function)

PCAs <- lapply(cell_type_DEG_output, function(x)(return(x$PCA)))
Volcanos <- lapply(cell_type_DEG_output, function(x)(return(x$Volcano)))

pdf("plots/DEG/PCA_compiled.pdf", width = 14, height = 10)
ggarrange(plotlist =  PCAs, common.legend = T)
dev.off()


pdf("plots/DEG/Volcanos_compiled.pdf", width = 15, height = 8)
ggarrange(plotlist =Volcanos, common.legend = T)
dev.off()

save(agg_cell_types, cell_type_DEG_output, de.results_act, file = "data/DEG.RData")


#### CHECK CELLTYPE ABUNDANCE DIFFERENCES ----
sce$broad_type <- "Somatic"
sce$broad_type[sce$new_clusters %in% c(
  "Spermatids", "Germ", "GSC and Spermatagonia")] <- "Germ"

cell_type_abundances_summary <- data.frame(type = sce$broad_type, sample = sce$sample, treatment = sce$treatment) %>% 
  dplyr::count(treatment, type, sample) %>% 
  spread(type,n)

glmer(cbind(cell_type_abundances_summary$Germ, 
            cell_type_abundances_summary$Germ + cell_type_abundances_summary$Somatic) ~ 
        cell_type_abundances_summary$treatment + (1|sample), data = cell_type_abundances_summary, family = "binomial") %>% 
  summary()


chisq_list <- list()
residuals <- list()
st <- cell_type_abundances_summary$sample[cell_type_abundances_summary$treatment == 'st']
sr <- cell_type_abundances_summary$sample[cell_type_abundances_summary$treatment == 'sr']

combos <- expand.grid(st, sr)
for (s in 1:nrow(combos)){
  tmp <- filter(cell_type_abundances_summary, sample %in% unlist(c(combos[s,1:2])))
  rownames(tmp) <- tmp$treatment
  cs <- chisq.test(tmp[,c(3,4)])
  chisq_list[[s]] <- cs
  residuals[[s]] <- cs$residuals
}

plot_res <- do.call(rbind, residuals) %>% as.data.frame()
plot_res$treatment <- rep(c("sr", "st"))

plot_res %>% 
  ggplot(aes(x = Somatic, colour = treatment)) + geom_density()


########## DELETE -------------
plot_data <- cell_type_abundances_summary %>% 
  mutate(Total = Germ + Somatic,
         Germ = Germ/Total, 
         Somatic = Somatic/Total) %>% 
  pivot_longer(c("Germ", "Somatic"), names_to = "Tissue",
               values_to = "Proportion of cell type")
plot_data[plot_data == 'sr'] <- "Drive"
plot_data[plot_data == 'st'] <- "Standard"
plot_data %>% 
  ggplot(aes(y = `Proportion of cell type`, x = treatment, fill = Tissue)) + geom_boxplot() + 
  theme_classic() 
ggsave("plots/sex_unfolded_figures/comp_cell_number.pdf", height = 4, width = 5)

