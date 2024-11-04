#rm(list = ls())
library(tidyverse)
library(ggrepel)
library(Seurat)
load("data/RData/DEG_DC.RData")
load("data/RData/seurat_final.RData")

ortholog_table <- read.csv("outdata/orthologs_Jan24.csv")
ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

MSLs1 <- c("PB.1152", "PB.9420")
MSLs2 <- c("gene-2030", "PB.520", "STRG.10819")
MSLs <- c(MSLs1, MSLs2)
MSLs <- c(MSLs1)


MSL_exp <- XA %>% 
  filter(treatment == "ST") %>% 
  filter(genes %in% MSLs)

med_msl <- medians %>% 
  merge(MSL_exp, by = c('celltype', 'treatment')) %>% 
  merge(ortholog_table[,c("REF_GENE_NAME", "consensus_gene")], by.x = 'genes', by.y = 'REF_GENE_NAME') %>% 
  filter(treatment == "ST") %>% unique()

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7")

p1 <- med_msl %>% 
  filter(treatment == "ST") %>% 
  ggplot(aes(x = median, y = logcpm, colour = celltype)) + geom_point() + 
  facet_wrap(~consensus_gene) + 
  xlab(median~(log[2]~X:A~expression~ratio)) + 
  ylab("log(CPM)") + 
  theme(legend.position = 'none') + 
  scale_colour_manual(values= cbPalette) + 
  geom_label_repel(data = med_msl, aes(x = median, y = logcpm, label = celltype),
                   box.padding   = 0.1, 
                   point.padding = 0.1, force = 10,
                   segment.color = 'grey50',show.legend = F, size = 2.8)

p2 <- DotPlot(seurat_final, features = MSLs,  assay = "RNA") + 
  labs(y = '', x = '') + 
  scale_x_discrete(labels = as.vector(c("msl-1", "msl-3", "mof", "Clamp", "mle"))) + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="right", 
        legend.box = "horizontol", 
        axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.ticks = element_line(color = "black")) + 
  coord_flip() + 
  guides(color = guide_colorbar(title = 'Scaled expression'), 
         size = guide_legend(title = '% of cells expressed'))

sup_msl <- cowplot::plot_grid(p1, p2, ncol = 1, axis = 'h', align = 'b',
                   labels = c("(a)", "(b)"), rel_heights = c(3,3))
sup_msl
sup_msl <- p2

ggsave("plots/MSL_expression.pdf", sup_msl, width = 6, height = 4)
system("open plots/MSL_expression.pdf")
