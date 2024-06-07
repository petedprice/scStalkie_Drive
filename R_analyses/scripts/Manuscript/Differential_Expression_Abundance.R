## based off http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#putting-it-all-together
# load libraries and functions
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


load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_cleaned_ss.RData")
Muscle <- c(17)
Spermatocytes <- c(10,11,13)
Spermatids <- c(3,4)
`GSC/Spermatogonia` <- c(0,2,7)
Pre_meiotic_cyst <- c(5,12,14,15)
Post_meiotic_cyst <- c(1,6,8,9,16)

### Plotting key markers against numbered cell clusters for assignment ---

seurat_integrated_ss@meta.data$celltype <- 'NA'
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Muscle] <- "Muscle"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatocytes] <- "Spermatocytes"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Spermatids] <- "Spermatids"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% `GSC/Spermatogonia`] <- "GSC & Spermatogonia"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Pre_meiotic_cyst] <- "Pre-meiotic \ncyst"
seurat_integrated_ss$celltype[seurat_integrated_ss$integrated_snn_res.0.4 %in% Post_meiotic_cyst] <- "Post-meiotic \ncyst"
seurat_integrated_ss$treatment <- "SR"
seurat_integrated_ss$treatment[grep("st", seurat_integrated_ss$sample)] <- "ST"
Idents(seurat_integrated_ss) <- seurat_integrated_ss$celltype

DefaultAssay(seurat_integrated_ss) <- "integrated"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]


##### BIOCONDUCTOR EDGER2 -----
## PLOT FUNCTION ---- 
make_pca_function <- function(x, ortholog_table){
  md <- metadata(x)
  mds_data <- md$y %>% 
    edgeR::cpm(log = TRUE) %>% 
    plotMDS()
  
  PCA <- data.frame(x = mds_data$x, 
                    y = mds_data$y, 
                    names = md$y$samples$seq_folder,
                    treatment = md$fit$samples$treatment) %>% 
    ggplot(aes(x = x, y = y, fill = names, label = names, colour = treatment))  +
    geom_text(hjust="inward", vjust="inward") +
    labs(x = "PC2", y = "PC1", title = md$y$samples$celltype[1]) 
  pdf(paste("plots/PCA_", md$y$samples$celltype[1], ".pdf", sep = ""))
  print(PCA) 
  dev.off()
  x <- x[is.na(x$logFC) == FALSE,] %>%
    as.data.frame()
  
  
  x$geneID <- gsub("-", "_", rownames(x))
  #assign x$gene using ortholog_table consensus_gene
  x <- merge(x, ortholog_table, by.x = 'geneID', by.y = 'REF_GENE_NAME')   
  top_p <- x$FDR[order(x$FDR)][5]
  top_p2 <- x$FDR[order(x$FDR)][5]
  x$celltype = md$y$samples$celltype[1]
  Volcano <- x %>% 
    ggplot(aes(x = logFC, y = -log10(FDR), label = consensus_gene,
               colour = (abs(logFC) > 1 &  FDR < 0.05))) + 
    geom_point() + geom_text_repel(data=subset(x, FDR < top_p & abs(logFC) > 1),
                                   aes(logFC,-log10(FDR),label=consensus_gene), 
                                   hjust="inward", vjust="inward") +
    labs(y = "-log10 Pvalue (FDR adjusted)", x = "log2 fold-change", 
         title = paste(md$y$samples$celltype[1], "(positive is ST bias)"))+ 
    theme(legend.position = "none")  + 
    geom_hline(yintercept = -log10(0.05), color = 'grey', linetype = 'dashed') + 
    geom_vline(xintercept = c(-1,1), color = 'grey', linetype = 'dashed')
  pdf(paste("plots/DEG_", md$y$samples$celltype[1], ".pdf", sep = ""))
  print(Volcano) 
  dev.off()
  ST_bias_genes <- filter(x, logFC > 1, FDR < 0.05)
  SR_bias_genes <- filter(x, logFC < -1, FDR < 0.05)
  output <- list(PCA,Volcano, ST_bias_genes, SR_bias_genes, x)
  names(output) <- c("PCA", "Volcano", "ST_bias", "SR_bias", "all")
  return(output)
  
}

sce <- #subset(seurat_integrated, celltype != 'Unknown') %>% 
  seurat_integrated_ss %>% as.SingleCellExperiment(assay = "RNA")


#### BULK RNASEQ ACROSS WHOLE TISSUE ----
whole_tissue <- aggregateAcrossCells(sce, id=colData(sce)[,c("sample")])
de.results_wt <- pseudoBulkDGE(whole_tissue, 
                               label='bulk',
                               condition = whole_tissue$treatment,
                               design=~treatment,
                               coef = "treatmentST",
)
metadata(de.results_wt$bulk)$y$samples$celltype <- "bulk"

bulk_DEG_output <- make_pca_function(de.results_wt$bulk, ortholog_table)



### PSEUDOBULK INDIVIDUAL CELL TYPES ----
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]

y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))
keep <- filterByExpr(y, group=agg_cell_types$treatment)
y <- y[keep,]
mds_data <-edgeR::cpm(y, log = TRUE) %>% 
  plotMDS()
plot <- data.frame(x = mds_data$x, 
                   y = mds_data$y, 
                   names = y$samples$seq_folder,
                   treatment = y$samples$treatment, 
                   cell_type = y$samples$celltype) %>% 
  ggplot(aes(x = x, y = y, label = cell_type, colour = treatment))  +
  geom_text(hjust="inward", vjust="inward") +
  labs(x = "PC2", y = "PC1", title = "ind cell types")

pdf("plots/PCA_pseudobulk_ind_cell_types.pdf")
plot 
dev.off()

de.results_act <- pseudoBulkDGE(agg_cell_types, 
                                label=agg_cell_types$celltype,
                                condition = agg_cell_types$treatment,
                                design=~treatment,
                                coef = "treatmentST"
)

cell_type_DEG_output <- lapply(de.results_act, make_pca_function, ortholog_table = ortholog_table)

PCAs <- lapply(cell_type_DEG_output, function(x)(return(x$PCA)))
Volcanos <- lapply(cell_type_DEG_output, function(x)(return(x$Volcano)))



pdf("plots/PCA_compiled.pdf", width = 14, height = 10)
ggarrange(plotlist =  PCAs, common.legend = T)
dev.off()


pdf("plots/Volcanos_compiled.pdf", width = 20, height = 15)
ggarrange(plotlist = Volcanos, common.legend = T)
dev.off()


dif_exp_data <- lapply(cell_type_DEG_output, function(x)(return(x$all))) %>% 
  bind_rows()

dif_exp_data$Significant <- ifelse(dif_exp_data$logFC > 1 & dif_exp_data$FDR < 0.05, "ST bias",
                                   ifelse(dif_exp_data$logFC < -1 & dif_exp_data$FDR < 0.05, "SR bias", "NS"))

dif_exp_figure <- dif_exp_data %>% 
  filter(Significant != "NS") %>% 
  ggplot(aes(x = celltype, fill = Significant)) + geom_bar(position = 'dodge') + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + labs(y = "Number of DEGs", 
                                                               x = 'Cell Type')



csq_chr_func <- function(x){
  data <- x$all
  data$sig <- 'Unbiased'
  data$sig[data$PValue < 0.05 & data$logFC > 1] <- 'ST bias'
  data$sig[data$PValue < 0.05 & data$logFC < -1] <- 'SR bias'
  chsq_table <- data %>%
    group_by(sig, chr) %>% 
    summarise(n = n()) %>% #trnasformn variables into rows and columns
    spread(sig, n)
  chsq_results <- chisq.test(chsq_table[,c(2:4)])
  rsds <- chsq_results$residuals %>% as.data.frame()
  rsds$chr <- chsq_table$chr
  rsds$celltype <- data$celltype[[1]]
  
  rsds$pvalue <- chsq_results$p.value[1]
  return(rsds)
} 

chisq_list <- lapply(cell_type_DEG_output, csq_chr_func) %>% 
  bind_rows() %>% #collapse not_sig, SR_bias and ST_bias into single column 
  pivot_longer(c("Unbiased", "SR bias", "ST bias"), names_to = "Expression class", values_to = "chisq_resid")


DGE_chisq_plot <- ggplot(chisq_list, aes(x = chr, y = chisq_resid, fill = `Expression class`)) + 
  facet_wrap(~celltype) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_classic() + 
 scale_fill_brewer(palette = "Set2") + 
  labs(x = "Chromosome", y = "Residuals")

ggarrange(DGE_chisq_plot, dif_exp_figure, ncol = 1, heights = c(2, 1), common.legend = T, 
          labels = c("A", "B"), legend = "bottom")
ggsave("plots/DGE_chisq_plot.pdf", width = 15, height = 20, units = "cm")

###################################################
########## DELETE ---------

#### CHECK CELLTYPE ABUNDANCE DIFFERENCES ----
abundances <- table(sce$celltype, sce$sample)
abundances <- unclass(abundances) 
head(abundances)

extra.info <- colData(sce)[match(colnames(abundances), sce$sample),]
y.ab <- DGEList(abundances, samples=extra.info)
keep <- filterByExpr(y.ab, group=y.ab$samples$treatment)
y.ab <- y.ab[keep,]
summary(keep)


design <- model.matrix(~factor(treatment), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)


plotBCV(y.ab, cex=1)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)


-------------------------------------------



### OLD ABUNDANCE DIFFERENCES --------
seurat_integrated$broad_type <- "Somatic"
seurat_integrated$broad_type[seurat_integrated$celltype %in% c(
  "Spermatids", "Spermatocytes", "GSC/Spermatogonia")] <- "Germ"
#sce$broad_type <- sce$celltype
seurat_integrated_ss <- subset(seurat_integrated, celltype %in% c("Cyst", "GSC/Spermatogonia", "Spermatids", "Spermatocytes"))
DimPlot(seurat_integrated_ss, group.by = "broad_type", label = T)
cell_type_abundances_summary <- data.frame(type = seurat_integrated_ss$broad_type, 
                                           sample = seurat_integrated_ss$sample, 
                                           treatment = seurat_integrated_ss$treatment) %>% 
  dplyr::count(treatment, type, sample) %>% 
  spread(type,n)
cell_type_abundances_summary$prop <- cell_type_abundances_summary$Germ / 
  (cell_type_abundances_summary$Germ + cell_type_abundances_summary$Somatic)


glmer(cbind(cell_type_abundances_summary$Germ, 
            cell_type_abundances_summary$Germ + cell_type_abundances_summary$Somatic) ~ 
        cell_type_abundances_summary$treatment + (1|sample), data = cell_type_abundances_summary, family = "binomial") %>% 
  summary()


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
plot_data %>% 
  ggplot(aes(y = `Proportion of cell type`, x = treatment, fill = Tissue)) + geom_boxplot() + 
  theme_classic() 
ggsave("plots/sex_unfolded_figures/comp_cell_number.pdf", height = 4, width = 5)


sum_data <- seurat_integrated_ss@meta.data %>% 
  group_by(celltype, treatment, sample) %>% 
  summarise(n = n()) %>% 
  spread(celltype,n)

sum_data$total <- rowSums(sum_data[,c(3:6)])
sum_data[,3:6] <- sum_data[,3:6]/sum_data$total

sum_data %>% 
  pivot_longer(c("Cyst", "GSC/Spermatogonia", "Spermatocytes", "Spermatids"), names_to = "celltype", values_to = "prop") %>%
  ggplot(aes(x = treatment, y = prop, fill = celltype)) + geom_boxplot() +
  theme_classic()
x <- cell_type_DEG_output$GSC_and_Spermatogonia
csq_chr_func <- function(x){
  data <- x$all
  data$sig <- 'Unbiased'
  data$sig[data$PValue < 0.05 & data$logFC > 1] <- 'ST_bias'
  data$sig[data$PValue < 0.05 & data$logFC < -1] <- 'SR_bias'
  chsq_table <- data %>%
    group_by(sig, chr) %>% 
    summarise(n = n()) %>% #trnasformn variables into rows and columns
    spread(sig, n)
  chsq_results <- chisq.test(chsq_table[,c(2:4)])
  rsds <- chsq_results$residuals %>% as.data.frame()
  rsds$chr <- chsq_table$chr
  rsds$celltype <- data$celltype[[1]]
  
  rsds$pvalue <- chsq_results$p.value[1]
  return(rsds)
} 

chisq_list <- lapply(cell_type_DEG_output, csq_chr_func) %>% 
  bind_rows() %>% #collapse not_sig, SR_bias and ST_bias into single column 
  pivot_longer(c("Unbiased", "SR_bias", "ST_bias"), names_to = "Expression class", values_to = "chisq_resid")


ggplot(chisq_list, aes(x = chr, y = chisq_resid, fill = `Expression class`)) + 
  facet_wrap(~celltype) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_classic() + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Chromosome", y = "Residuals")



