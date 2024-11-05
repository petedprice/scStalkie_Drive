make_plots_function <- function(x, plotpath, name = "metadata_QC"){
  print('yep')
  plots <- list()
  plots[[1]] <- ggplot(x, aes(x=sample, fill=sample)) + geom_bar() +theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  plots[[2]] <- ggplot(x, aes(color=sample, x=nUMI, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  plots[[3]] <- ggplot(x, aes(color=sample, x=nGene, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  plots[[4]] <- ggplot(x, aes(x=sample, y=log10(nGene), fill=sample)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells vs NGenes")
  plots[[5]] <- ggplot(x, aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~sample)
  plots[[6]] <- ggplot(x, aes(color=sample, x=mitoRatio, fill=sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  plots[[7]] <- ggplot(x, aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  comp_plots <- ggarrange(plotlist = plots)
  comp_plots <- annotate_figure(comp_plots, top = text_grob(name,
                                                            color = "Black", face = "bold", size = 14))
  ggsave(filename = paste(plotpath, "/", name, ".pdf", sep = ""), comp_plots, width = 30, height = 20)

}



marker_gene <- function(seurat_object){
  
}

