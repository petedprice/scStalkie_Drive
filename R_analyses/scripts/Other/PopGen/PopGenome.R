library(PopGenome)
library(tidyverse)


CHRX <- read.csv("indata/PopGen/Chr_X_PopGenome_Stats/default_Chr_X_PopGenome_Stats.csv", header = T)
CHRX$chr <- "Chr_X"
#CHR1 <- read.csv("indata/PopGen/PopGenome_Stats/Chr_1_PopGenome_Stats/default_Chr_1_PopGenome_Stats.csv", header = T)
#CHR1$chr <- "Chr_1"
#CHR2 <- read.csv("indata/PopGen/PopGenome_Stats/Chr_2_PopGenome_Stats/default_Chr_2_PopGenome_Stats.csv", header = T)
#CHR2$chr <- "Chr_2"

#all_chroms <- rbind(CHRX, CHR1, CHR2)
all_chroms <- CHRX
all_chroms$Pi_alls <- all_chroms$Pi_all/all_chroms$n.sites_all
all_chroms$STs_pis <- all_chroms$STs_pi/all_chroms$ST_n.segregating.sites
all_chroms$SRs_pis <- all_chroms$SRs_pi/all_chroms$SR_n.segregating.sites

all_chroms %>% head()
all_chroms$Pos <- all_chroms$Pos/10000000
all_chroms %>% 
  filter(ST_n.segregating.sites > 9) %>% 
  filter(SR_n.segregating.sites > 9) %>%
  ggplot(aes(x = Pos, y = FST)) +geom_line() + 
  xlim(0,7)
 # geom_point(alpha = 0.1) + 
#  facet_wrap(~chr)
all_chroms %>% 
  ggplot(aes(x = Pos, y = STs_pis)) + geom_point()

vcf <- read.table("indata/PopGen/Chr_1_filt_fin.vcf.gz", sep = '\t', header = T, comment.char = '#')

ggplot(vcf, aes(x = X56)) + geom_density()
all_chroms %>% 
  ggplot(aes(x = chr, y = log(SR_n.segregating.sites))) + 
  geom_boxplot()

all_chroms %>%
  ggplot(aes(x = chr, y = log(ST_n.segregating.sites))) +
  geom_boxplot()

head(vcf)
vcf %>% 
  filter(X56 > 5150000 & X56 < 5160000)  %>% View()#%>% 
#  ggplot(aes(x = X56)) + geom_density()





CHRX <- read.csv("indata/PopGen/Chr_X_PopGenome_Stats/default_Chr_X_genelevel_PopGenome_Stats.csv", header = T)
CHRX$chr <- "Chr_X"

CHRX %>% 
  ggplot(aes(x = fisher.P.value, y = neutrality.index)) + geom_point()




genes <- unique(ortholog_table$REF_GENE_NAME)
CHRX_genes <- vector(length = ncol(CHRX) + 1)
for (gene in genes){
  print(gene)
  tmp <- filter(CHRX, as.numeric(Pos) >= ortholog_table$start[ortholog_table$REF_GENE_NAME == gene] & 
                  as.numeric(Pos) <= ortholog_table$end[ortholog_table$REF_GENE_NAME == gene])
  if (nrow(tmp) >0){
    tmp$gene <- gene
    CHRX_genes <- rbind(CHRX_genes, tmp)
  }
}

CHRX_genes <- CHRX_genes[-1,]

CHRX_genes_DGE <- CHRX_genes %>% 
  merge(dif_exp_data, by.x = 'gene', by.y = 'genes')


sumr <- CHRX_genes_DGE %>% 
  group_by(gene, FDR, logFC, celltype, Significant) %>% 
  summarise(FST = mean(FST), 
            SR_D = mean(SR_Tajima.D), 
            ST_D = mean(ST_Tajima.D), 
            SR_seg = mean(SR_n.segregating.sites), 
            ST_seg = mean(ST_n.segregating.sites), 
            SR_pi = mean(SRs_pi), 
            ST_pi = mean(STs_pi))


sumr %>% 
  ggplot(aes(x = SR_pi, y = log(abs(logFC)))) + geom_point() +
  facet_wrap(~celltype) +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson')

sumr %>% 
  ggplot(aes(x = FST, y = log(abs(logFC)))) + geom_point() + 
  facet_wrap(~celltype) + 
  geom_smooth(method = 'lm') + 
  ggpubr::stat_cor(method = 'pearson')


sumr %>% 
  ggplot(aes(x = SR_D, y = log(abs(logFC)))) + geom_point() +
  facet_wrap(~celltype) + 
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson')

sumr %>% 
  ggplot(aes(x = ST_D, y = log(abs(logFC)))) + geom_point() +
  facet_wrap(~celltype) + 
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson')


sumr %>% 
  filter(Significant != "Unbiased") %>%
  ggplot(aes(x = Significant, y = SR_D)) + geom_boxplot() +
  facet_wrap(~celltype)


sumr %>% 
  filter(Significant != "Unbiased") %>%
  ggplot(aes(x = Significant)) + geom_bar() + 
  facet_wrap(~celltype)

