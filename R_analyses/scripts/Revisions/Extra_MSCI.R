
load("data/RData/DEG_DC.RData")


XA %>% 
  filter(Chr != "Mito") %>% 
  ggplot(aes(x = treatment, y = logcpm, fill = Chr)) + 
  geom_boxplot() + 
  facet_wrap(~celltype)
