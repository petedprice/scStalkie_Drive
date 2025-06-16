rm(list = ls())
library(tidyverse)
library(ggrepel)
library(Seurat)
load("data/RData/DEG_DC.RData")
load("data/RData/seurat_final.RData")

ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")

ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]
ortholog_table$REF_GENE_NAME <- gsub("_", "-", ortholog_table$REF_GENE_NAME)
ortholog_table$consensus_gene <- gsub("_", "-", ortholog_table$consensus_gene)

MSLs1p <- c("PB.1152", "PB.9420", "gene-2030", "PB.520", "STRG.10819", "gene-5827")
MSL_df1 <- ortholog_table %>% 
  dplyr::select(REF_GENE_NAME, consensus_gene) %>% 
  filter(REF_GENE_NAME %in% MSLs1p) %>% 
  unique()
MSLs1v <- c("PB.1657", "msl-2")
MSL_df <- rbind(MSL_df1, MSLs1v)
MSL_df <- MSL_df[order(MSL_df$consensus_gene),]

ortholog_table %>% 
  filter(REF_GENE_NAME %in% MSL_df$REF_GENE_NAME) %>% 
  dplyr::select(REF_GENE_NAME, OMA_DMEL, OF_DMEL, consensus_gene) %>% 
  unique() %>% 
  rename(Tdal_gene = REF_GENE_NAME, 
         Orthofinder_ortholog = OF_DMEL, 
         OMA_ortholog = OMA_DMEL) %>% 
  write.csv("data/DCC.csv")


sup_msl <- DotPlot(seurat_final, features = MSL_df$REF_GENE_NAME,  assay = "RNA") + 
  labs(y = '', x = '') + 
  scale_x_discrete(labels = as.vector(c(MSL_df$consensus_gene))) + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="right", 
        legend.box = "horizontol", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) + 
  coord_flip() + 
  guides(color = guide_colorbar(title = 'Scaled expression'), 
         size = guide_legend(title = '% of cells expressed'))


ggsave("plots/manuscript_plots/S5.tiff", sup_msl, width = 6, height = 4)
system("open plots/manuscript_plots/S5.tiff")
