######################### FIGURE 1 (and old 2 incorporated) ######################### 

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7")

expected_X_exp <- filter(ortholog_table, chr == "Chr_X") %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
expected_A_exp <- filter(ortholog_table, chr %in% c("Chr_1", "Chr_2")) %>%
  dplyr::select(REF_GENE_NAME) %>% 
  unique() %>% nrow()
expected_expression = expected_X_exp/(expected_A_exp + expected_X_exp)

dps_all_figure <- DotPlot(seurat_final, features = names(mfms), assay = "RNA")+coord_flip() +
  scale_x_discrete(labels = as.vector((mfms))) + 
  labs(y = '', x = "Marker gene expression") + 
  theme_classic() + scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="bottom", legend.box = "vertical", 
        axis.text.x = element_text(angle = 45, hjust=1))
    



nfeatures_figure <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA_autosomes, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + scale_fill_manual(values= cbPalette) + 
  labs(x = "", y = "N\u00b0 detected autosomal genes") + 
  theme(legend.position="none", 
        axis.line.x = element_blank()) + rm_axis


UMAP <- DimPlot(seurat_final, cols = cbPalette, group.by = 'celltype', label = T) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(legend.position="none")
UMAP2 <- UMAP[[1]]$data %>% 
  ggplot(aes(x = umap_1, y = umap_2, colour = celltype)) + 
  geom_point(stroke = NA, size = .7) + scale_color_manual(values= cbPalette) + 
  guides(colour = guide_legend(override.aes = list(size=3), title = "Cell type")) + 
  theme_classic() + labs(x = "UMAP 1", y = "UMAP 2", legend = 'Cell type') +
  theme(legend.position = 'none') +  
  shadowtext::geom_shadowtext(data = UMAP[[1]]$layers[[2]]$data, aes(label = celltype), colour = 'black', 
                              size = 4, nudge_y = -0.5, nudge_x = 0, 
                              bg.color = 'white', bg.r = 0.1, bg.size = 0.5, 
                              fontface = "bold") + 
  geom_point(data = UMAP[[1]]$layers[[2]]$data, aes(x = umap_1, y = umap_2+0.1), colour = 'black',
             stroke = 1, size = 2, shape = 21, fill = 'white')

# DC Plot #
check_normalities <- 
  XA %>% mutate(ratio = XA) %>% 
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  group_by(celltype, treatment) %>% 
  summarise(p_value = shapiro.test(ratio)$p.value)

#Ratios not normally distributed 

# Perform t-test for each group and store results in a data frame
ST_wilcox_results <- XA %>%
  filter(Chr == "X") %>% 
  filter(logcpm >= 1) %>% 
  filter(treatment == "ST") %>% 
  group_by(celltype, treatment) %>%
  summarise(
    p_value = wilcox.test(XA, mu = 0, alternative = 'two.sided')$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.00001 ~ "***",
      p_value < 0.001 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ " "
    )
  ) %>% 
  mutate(treatment = "ST")


DC_plot <- XA %>% 
  filter(Chr == "X" & logcpm >= 1 & treatment == "ST") %>% 
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Spermatids"))) %>% 
  mutate(Chromosome = Chr) %>% 
  ggplot(aes(x = celltype, y = (XA), fill = treatment)) + geom_boxplot(outlier.alpha = 0.2) + 
  labs(x = "", y = "log(CPM)", fill = "Chromosome") + #add sig values between Chrs
  geom_hline(yintercept = c(0,-1), linetype = 'dashed', colour = 'black', linewidth = 0.7) + 
  
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  scale_y_continuous(breaks=c(0.0, 5.0, 10.0, 15.0, 20.0)) + 
  guides(fill = 'none') +
  ylab(expression(log[2]~X:A~ratio)) + 
  geom_text(data = ST_wilcox_results, aes(x = celltype, y = 9, label = significance),
            vjust = -0.5, size = 5) + 
  ylim(-5,10) + rm_axis + 
  theme(axis.line.x = element_blank())

# XCI Plot #

XCI_wilcox_results <- seurat_final@meta.data %>%
  filter(treatment == "ST") %>% 
  group_by(celltype, treatment) %>%
  summarise(
    p_value = wilcox.test(Xprop, mu = expected_expression, alternative = 'two.sided')$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.00001 ~ "***",
      p_value < 0.001 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ " "
    )
  )


XCI_plot <- seurat_final@meta.data %>% 
  filter(treatment == "ST") %>% 
  ggplot(aes(x = celltype, y = Xprop, fill = treatment)) + geom_violin() +
  theme_classic() + scale_fill_brewer(palette = 'Set2') + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill = 'none') +
  labs(y = "N\u00b0 expressed X genes / \n N\u00b0 expressed genes", x = "")

# Combine all plots #

A <- ggplotGrob(UMAP2)
B <- ggplotGrob(dps_all_figure)
C <- ggplotGrob(nfeatures_figure)
D <- ggplotGrob(DC_plot)
E <- ggplotGrob(XCI_plot)

AB <- align_plots(A,B, align = 'v', axis = 'l')
DE <- align_plots(D,E, align = 'v', axis = 'l')
BE <- plot_grid(AB[[2]],DE[[2]], align = 'h', axis = 'bt', labels = c("(b)", "(e)"), rel_widths = c(2,1))
CD <- plot_grid(C,DE[[1]], align = 'v', axis = 'l', nrow = 2, labels = c("(c)", "(d)"))
ACD <- plot_grid(AB[[1]], CD, align = 'h', axis = 'b', ncol = 2, labels = c("(a)", ""), rel_widths = c(2,1))
Fig1 <- plot_grid(ACD, BE, align = 'hv', axis = 'l', nrow = 2, rel_heights = c(3,2))

ggsave("plots/FIG1.pdf", Fig1, width = 13, height = 13)
system("open plots/FIG1.pdf")

### FIGURE S1 ###

S1a <- DimPlot(seurat_final, group.by = 'Phase') + 
  labs(title = "Cell-cycle Phase") + 
  labs(x = "UMAP 1", y = "UMAP 2")

S1b <- DotPlot(seurat_final, features = names(afms), assay = "RNA")+
  scale_x_discrete(labels = as.vector((afms))) + coord_flip() +
  labs(y = '', x = "Marker genes expression") + 
  theme_classic() + 
  scale_colour_gradientn(colours = colorspace::diverge_hcl(8)) + 
  theme(legend.position="right", 
        legend.box = "horizontol",
        axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title=element_text(size=13))

S1c <- seurat_final@meta.data %>% 
  ggplot(aes(x = celltype, y = nFeature_RNA, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.2) + 
  theme_classic() + 
  scale_fill_manual(values= cbPalette) + 
  labs(x = "", y = "N\u00b0 detected genes") + 
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title=element_text(size=13))


FigS1 <- plot_grid(S1a, plot_grid(S1b + rm_axis, S1c, nrow = 2, 
                                  align = 'v', labels = c("B", "C"), 
                                  rel_heights = c(2.2,2)), 
                   nrow = 1, align = 'h', rel_widths = c(1, 0.7), axis = 't', 
                   labels =c("A", "")) 
ggsave("plots/FIGS1_MarkersUMAPFeatures.pdf", FigS1, height = 9, width =20)
system("open plots/FIGS1_MarkersUMAPFeatures.pdf")

#####################################################


