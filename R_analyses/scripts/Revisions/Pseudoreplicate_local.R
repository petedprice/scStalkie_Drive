library(tidyverse)
library(ggpubr)

rm(list = ls())

summarys <- lapply(list.files("data/pseudoreplicate/sumarised/", full.names = T, pattern = ".csv.gz"), read.csv) %>% 
  bind_rows()


#Summarie for both directions of DGE
summarys <- summarys %>% 
  mutate(sig = ifelse(sig != "non", "sig", "non")) %>% 
  group_by(chr, celltype, treatment, rep, sig) %>%
  reframe(prop = sum(prop), 
            n = sum(n)) %>% 
  group_by(chr, celltype, treatment, rep) %>% 
  mutate(celltype = ifelse(celltype == "GSCandSpermatogonia", "GSC/Spermatogonia", celltype),
         celltype = gsub("_", " ", celltype),
         celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")),
         chr = ifelse(chr != "Chr_X", "Autosomes", "X chromosome"))   



#Get p values for cmoparing groups for each celltype
p_val_across_chr <- summarys %>% 
  mutate(sig = ifelse(sig == "non", "non", "sig")) %>% 
  group_by(celltype, treatment, rep, sig) %>% 
  summarise(n = sum(n)) %>% 
  ungroup(sig) %>% 
  mutate(prop = n/sum(n)) %>% 
  filter(sig == 'sig') %>% 
  group_by(celltype) %>% 
  summarise(significant = wilcox.test(prop[treatment == "STvsSR"], 
                                    prop[treatment == 'STvsST'], alternative = 'two.sided')$p.value) %>% 
  mutate(group1 = "STvsSR", 
         group2 = "STvsST", 
         y.position = 30) %>% 
  mutate(sig = ifelse(significant > 0.05, "ns", ifelse(
    significant > 0.01, "*", ifelse(
      significant > 0.001, "**", "***")
    )
  ))

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "yellow3", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = cbPalette, 
                                                                  alpha = 0.5))

plot <- summarys %>% 
  mutate(prop = prop * 100) %>% 
  filter(sig != "non") %>% 
  ggplot(aes(x = chr, fill = treatment, y = prop)) + geom_boxplot() + 
  stat_compare_means(method = 'wilcox', aes(label = after_stat(p.signif)), label.y = 50) + 
  scale_fill_grey(start = 0.4, end = 0.7) + 
  theme_classic() + 
  labs(fill = "pseudo-replicate comparison", x = "", y = "% of genes differentially expressed") + 
  ylim(0,65) + theme(legend.position = 'top') + 
  ggh4x::facet_wrap2(~ celltype, strip = strip, ncol = 2)



ggsave("plots/manuscript_plots/S9.tiff", plot, height = 8, width = 6)
system("open plots/manuscript_plots/S9.tiff")



  
   















### DELETE OLD 

plot <- summarys %>% 
  mutate(prop = prop * 100) %>% 
  filter(sig != "non") %>% 
  ggboxplot(., x = 'treatment', y = 'prop', fill = 'chr') + 
  stat_pvalue_manual(p_val_across_chr, label = 'sig', y.position = 60) + 
  stat_compare_means(aes(group = chr, label = after_stat(p.signif)), method = 'wilcox.test', label.y  = 48) + 
  facet_wrap(~celltype, ncol = 2) + 
  scale_fill_grey(start = 0.4, end = 0.7) + 
  theme_classic() + 
  labs(fill = "", x = "pseudo-replicate comparison", y = "% of genes differentially expressed") + 
  ylim(0,65)



summarys %>% filter(sig != "non") %>% 
  ggplot(aes(x = treatment, y = prop * 100, fill = chr)) + 
  geom_boxplot(outlier.alpha = 0.011, position = position_dodge(0.75)) +
  scale_fill_grey(start = 0.4, end = 0.7) + 
  theme_classic() + 
  facet_wrap(~celltype) + 
  labs(fill = "", y = "% of genes differentially expressed") + 
  stat_compare_means(method = "wilcox.test", 
                     aes(
                       #label = after_stat(p.signif)
                     ), 
                     label.y = 40) 
#  ylim(0,45)# + 
#  stat_pvalue_manual(p_val_across_chr, label = 'significant', aes(group =))

#stat_pvalue_manual(p_values[p_values$chr == "Chr_X"&p_values$p < 0.05,], label = 'sig_fdr') + 

