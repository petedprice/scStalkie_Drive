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


load("data/RData/seurat_final.RData")
DefaultAssay(seurat_final) <- "RNA"
ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

Xgenes <- filter(ortholog_table, chr == "Chr_X")$REF_GENE_NAME %>% 
  intersect(rownames(seurat_final))
Agenes <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2"))$REF_GENE_NAME %>% 
  intersect(rownames(seurat_final))
sce <- seurat_final %>% as.SingleCellExperiment(assay = "RNA")

################################################################################



####################### PLOT PCAS and DGE FUNCTION #############################
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

################################################################################



############################ BULK RNASEQ ACROSS WHOLE TISSUE ###################
whole_tissue <- aggregateAcrossCells(sce, id=colData(sce)[,c("sample")])
de.results_wt <- pseudoBulkDGE(whole_tissue, 
                               label='bulk',
                               condition = whole_tissue$treatment,
                               design=~treatment,
                               coef = "treatmentST",
)
metadata(de.results_wt$bulk)$y$samples$celltype <- "bulk"

bulk_DEG_output <- make_pca_function(de.results_wt$bulk, ortholog_table)

################################################################################



######################## PSEUDOBULK INDIVIDUAL CELL TYPES ######################
agg_cell_types <- aggregateAcrossCells(sce, id=colData(sce)[,c("celltype", "sample")])
agg_cell_types <- agg_cell_types[,agg_cell_types$ncells >= 10]
y <- DGEList(counts(agg_cell_types), samples=colData(agg_cell_types))

design <- model.matrix(~agg_cell_types$treatment+agg_cell_types$celltype)
keep <- filterByExpr(y, design=design)
y <- y[keep,]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
cpm <- edgeR::cpm(y, log=T) %>% as.data.frame()

mds_data <-cpm %>% 
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


################################################################################




########################## DOSAGE COMPENSATION PLOTS ---------------------------

raw_counts <- y$counts %>% as.data.frame()
cpm$genes <- rownames(cpm)
raw_counts$genes <- rownames(raw_counts)

cols <- c("treatment", "celltype", "sample")
tidy_cpm <- cpm %>% as.data.frame %>% 
  pivot_longer(cols = colnames(cpm)[1:(length(colnames(cpm))-1)], 
               names_to = "sample", values_to = "logcpm") %>% 
  merge(y$samples[,cols], by.x = "sample", by.y = 0) %>% 
  rename(edgeRsample = sample,  sample = sample.y) %>% 
  mutate(Chr = ifelse(genes %in% Xgenes, "X", "A")) %>% 
  mutate(Chr = ifelse(!(genes %in% c(Xgenes, Agenes)), "Mito", Chr))

tidy_raw_counts <- raw_counts %>% as.data.frame %>% 
  pivot_longer(cols = colnames(raw_counts)[1:(length(colnames(raw_counts))-1)], 
               names_to = "sample", values_to = "raw_counts") %>% 
  merge(y$samples[,cols], by.x = "sample", by.y = 0) %>% 
  rename(edgeRsample = sample,  sample = sample.y) %>% 
  mutate(Chr = ifelse(genes %in% Xgenes, "X", "A")) %>% 
  mutate(Chr = ifelse(!(genes %in% c(Xgenes, Agenes)), "Mito", Chr))

check_dc_plot <- tidy_cpm %>% 
  filter(logcpm > 2) %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Pre-meiotic \ncyst", 
                                                "Post-meiotic \ncyst", "GSC/Spermatogonia", 
                                                "Primary Spermatocytes", "Secondary Spermatocytes", 
                                                "Spermatids"))) %>%
  filter(Chr != 'Mito') %>% 
  group_by(celltype, genes, Chr, treatment) %>% 
  summarise(logcpm = mean(logcpm)) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = logcpm, fill = Chromosome)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "Cell Type", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 20) + 
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  facet_wrap(~treatment)


ggsave("plots/dosage_compensation.pdf", check_dc_plot, width = 10, height = 6)

################################################################################



############################## CELLTYPE DGE ------------------------------------

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

################################################################################



######################## CHISQUARED ENRICHMENT FIGURE ##########################

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

################################################################################




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



