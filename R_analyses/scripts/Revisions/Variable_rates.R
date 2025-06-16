library(tidyverse)
library(ape)
library(stringr)
library(plyranges) 
library(rtracklayer)
library(ggpubr)
rm(list = ls())

load("data/RData/seurat_final.RData")
load("data/RData/DEG_DC.RData")
load("data/trajectory/sce_GAMed_8k.RData")
load("data/RData//traj_analyis.RData")
seurat_final <- JoinLayers(seurat_final)
Idents(seurat_final) <- "celltype"
levels(seurat_final) <- c("Muscle", "Early cyst", "Late cyst", 
                          "GSC/Spermatogonia", "Primary spermatocytes", 
                          "Secondary spermatocytes", "Early spermatids", "Late spermatids")
seurat_final@meta.data <- seurat_final@meta.data %>%  
  mutate(celltype = factor(celltype, levels = c("Muscle", "Early cyst", 
                                                "Late cyst", "GSC/Spermatogonia", 
                                                "Primary spermatocytes", "Secondary spermatocytes", 
                                                "Early spermatids", "Late spermatids")))

DefaultAssay(seurat_final) <- "RNA"

ortholog_table <- read.table("outdata/orthologs_April25.tsv", sep = '\t', header = T, 
                             stringsAsFactors = F, quote = "", comment.char = "")


ortholog_table$consensus_gene[is.na(ortholog_table$consensus_gene)] = 
  ortholog_table$REF_GENE_NAME[is.na(ortholog_table$consensus_gene)]

species <- c("T_dal")
names(species) <- c("Teleopsis_dalmanni")


branch_folders <- list.files("data/var_rates/gene_models/", full.names = T,  
                             pattern = '_branch_models')

branch_files <- list.files(branch_folders, full.names = T, recursive = T)
x <- branch_files[1]
get_omega_func <- function(x){
  split1 <- str_split(x, "/")[[1]]
  sp <- str_split(split1[length(split1)], "_")[[1]][1:2] %>% paste0(collapse = "_")
  orthogroup <- str_split(split1[length(split1)], "_")[[1]][3]
  #sp_name <- species[sp]
  sp_name <- sp
  x_lines <- readLines(x)
  start <- which(x_lines == 'dN & dS for each branch')+4
  end <- grep("tree length for dN:", x_lines)-2
  indx <- c(start:end)
  
  data <- x_lines[indx] %>% 
    strsplit(., " ") %>% 
    sapply(., function(x) x[x != ""]) %>% t() %>% 
    as.data.frame() %>% 
    dplyr::rename(branch = V1, t = V2, N = V3, S = V4, `dNdS` = V5, 
                  dN = V6, dS = V7, `N*dN` = V8, `S*dS` = V9)
  data$species_node <- str_split(data$branch, "\\.", simplify = T)[,3]
  data$species_node_name <- "node"
  
  branch_check_indx <- grep('tree length =', x_lines)+2
  branch_check <- x_lines[c(branch_check_indx, branch_check_indx+2)] %>% 
    strsplit(., "")
  
  bc_numbers <- branch_check[[1]]
  bc_names <- branch_check[[2]]
  writeLines(bc_numbers, "tmp_bc_numbers.txt")
  writeLines(bc_names, "tmp_bc_names.txt")
  
  bc_numbers_tree <- read.tree("tmp_bc_numbers.txt")$tip.label
  bc_names_tree <- read.tree("tmp_bc_names.txt")$tip.label
  names(bc_names_tree) <- bc_numbers_tree
  #use this naming to change species_node_name 
  data$species_node_name <- bc_names_tree[data$species_node]
  
  species_dNdS <- data %>% 
    filter(species_node_name == sp_name) %>% 
    select(dNdS) %>% 
    as.numeric()
  
  ssdata <- data %>% 
    filter(as.numeric(dNdS) == as.numeric(species_dNdS))
  species_dS <- sum(as.numeric(ssdata$dS))
  species_dN <- sum(as.numeric(ssdata$dN))
  species_N <- sum(as.numeric(ssdata$N))
  species_S <- sum(as.numeric(ssdata$S))
  NdN <- sum(as.numeric(ssdata$N) * as.numeric(ssdata$dN))
  SdS <- sum(as.numeric(ssdata$S) * as.numeric(ssdata$dS))
  out_data <- data.frame(species = sp, species_sci = sp_name, dS = species_dS, 
                         dN = species_dN, NdN = NdN, SdS = SdS, N = species_N, S = species_S,
                         dNdS = species_dNdS, orthogroup = orthogroup)
}


omegas <- lapply(branch_files, get_omega_func) %>% 
  bind_rows()
rownames(omegas) <- NULL

get_name_func <- function(x){
  str_split(x, " ")[[1]][1:2] %>% gsub(">", "", .) %>% 
    gsub("gene=", "", .)
}  

gene_names <- readLines("data/var_rates/Teleopsis_dalmanni.cds_longest.fna") %>% 
  .[grepl(">", .)] %>% 
  lapply(., get_name_func) %>% 
  do.call(rbind, .) %>% as.data.frame() %>% 
  dplyr::rename(longest_transcript = V1, genes = V2)

dif_exp_sum <- dif_exp_data %>% 
  group_by(genes, chr, consensus_gene) %>% 
  summarise(ST_bias = sum(Significant == "ST-biased") > 0, 
            SR_bias = sum(Significant == "SR-biased") > 0, 
            logFC_DGE = mean(logFC)) %>% 
  mutate(DGE_bias = ifelse(ST_bias == T & SR_bias == T, "BOTH", 
                           ifelse(ST_bias == T, "ST", 
                                  ifelse(SR_bias == T, "SR", "UNBIAS")))) %>% 
  dplyr::select(-c("ST_bias", "SR_bias"))


orthogroups <-read.table("data/var_rates/N0.tsv", sep = "\t", header = T)
omega_genes <- merge(omegas, orthogroups, by.x = "orthogroup", by.y = 'OG', all.x = T, all.y = F) %>% 
  merge(gene_names, by.x = 'Teleopsis_dalmanni.protein_longest', by.y = 'longest_transcript') %>% 
  mutate(genes = gsub("_", "-", genes)) %>% 
   merge(dif_exp_sum, by = 'genes')
 


## APPROACH 1

output <- list()
for (class in unique(omega_genes$DGE_bias)){
  genes <- filter(omega_genes, DGE_bias == class)$genes
  for (i in 1:1000){
    gene_sample <- genes %>%sample(., 10, replace = T) %>% unique()
    genes_df <- omega_genes %>% 
      filter(genes %in% gene_sample) %>% 
      filter(!is.na(dNdS)) %>%
      summarise(S = sum(S, na.rm = T), 
                N = sum(N, na.rm = T),
                SdS = sum(SdS, na.rm = T),
                NdN = sum(NdN, na.rm = T),
                dS = sum(dS, na.rm = T), 
                dN = sum(dN, na.rm = T),
                dS2 = SdS/S,
                dN2 = NdN/N) %>% 
      mutate(dNdS = dN/dS, 
             dNdS2 = dN2/dS2) %>% 
      mutate(class = class)
    output[[paste0(i, class)]] <- genes_df
  }
}

output <- do.call(rbind, output) 

comparisons = list(c("SR", "ST"), c("ST", "UNBIAS"), c("SR", "UNBIAS"))


a <- output %>% 
  ggplot(aes(x = class, y = dNdS, fill = class)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  stat_compare_means(comparisons = comparisons, method = 'wilcox') + 
  theme_classic() + scale_fill_grey(start = 0.3, end = 0.7) + 
  theme(legend.position = 'none')




b <- omega_genes %>% 
  #filter(DGE_bias != "UNBIAS") %>% 
    filter(dS <=2) %>%
    filter(dNdS < 1) %>% 
  dplyr::select(genes, dNdS, dS, dN, consensus_gene, DGE_bias) %>%
  unique() %>% 
  ggplot(aes(x = DGE_bias, y = dNdS, fill = DGE_bias)) +
  geom_boxplot() + 
  stat_compare_means(comparisons = comparisons, method = 'wilcox') + 
  theme_classic() + scale_fill_grey(start = 0.3, end = 0.7) + 
  theme(legend.position = 'none')



plot <- ggarrange(a,b, labels = c("(A)", "(B)"), label.x = -0.03)
ggsave("plots/manuscript_plots/variable_rates.pdf", plot, height = 5, width = 10)
system('open plots/manuscript_plots/variable_rates.pdf')





############### APPROACH 2 
omega_genes <- omega_genes %>% 
  mutate(chr2 = ifelse(chr == "Chr_X", "Chr_X", "Chr_A"), 
         chr = ifelse(chr == "Chr_X", "Chr_X", "Chr_A"))

#Confidence intervals 
library(boot)

dNdS_func <- function(df, index){
  out <- df[index,] %>% 
    summarise(S = sum(S, na.rm = T), 
                   N = sum(N, na.rm = T),
                   SdS = sum(SdS, na.rm = T),
                   NdN = sum(NdN, na.rm = T),
                   dS = sum(dS, na.rm = T), 
                   dN = sum(dN, na.rm = T),
                   dS2 = SdS/S,
                   dN2 = NdN/N) %>% 
    mutate(dNdS = dN/dS, 
           dNdS2 = dN2/dS2)
  return(out$dNdS)
}

output <- list()


output <- list()
for (Chrom in unique(omega_genes$chr2)){
  for (class in unique(omega_genes$DGE_bias)){
    ss <- filter(omega_genes, chr2 == Chrom & DGE_bias == class)
    boot_dist <- boot(data = ss, statistic = dNdS_func, R = 1000)
    boot_dist$class <- class
    boot_dist$chr <- Chrom
    output[[paste0(class, Chrom)]] <- boot_dist
    
  }
}

get_quantiles_func <- function(boot_dist){
  CIs = quantile(boot_dist$t,  probs=c(0,.025,.25,.75, .975,1))
  out <- data.frame(chr = boot_dist$chr, 
                    class = boot_dist$class, 
                    CI0 = CIs[1],
                    CI2.5=CIs[2],
                    CI25=CIs[3],
                    CI75=CIs[4],
                    CI97.5=CIs[5],
                    CI1=CIs[6],
                    mean = boot_dist$t0, 
                    max = boot_dist$t0, 
                    min = boot_dist$t0)
  row.names(out) <- NULL
  return(out)
}

bs_data_func <- function(boot_dist){
  out <- data.frame(values = boot_dist$t, 
                    chr = boot_dist$chr, 
                    class = boot_dist$class)
  return(out)
}


boostrap_data <- lapply(output, bs_data_func) %>% 
  bind_rows()



quantiles <- lapply(output, get_quantiles_func) %>% 
  bind_rows()

quantiles %>% 
  ggplot(aes(x = class, colour = chr, y = mean)) +
  geom_point() +  
  geom_errorbar(aes(ymin = CI0, ymax = CI1), colour = 'black') + 
  geom_errorbar(aes(ymin = CI2.5, ymax = CI97.5), colour = 'grey') +
  facet_wrap(~chr) + 
  theme_classic()


lines <- data.frame(comb = c("SR_ST", "SR_UNBIAS", "ST_UNBIAS", "SR_ST", "SR_UNBIAS", "ST_UNBIAS"),
                    chr = c("Chr_A", 'Chr_A', "Chr_A", "Chr_X", "Chr_X", "Chr_X"),
                    true_dif = c(quantiles$mean[3]- quantiles$mean[2],
                            quantiles$mean[3]- quantiles$mean[1],
                            quantiles$mean[2]- quantiles$mean[1],
                            quantiles$mean[6]- quantiles$mean[5],
                            quantiles$mean[6]- quantiles$mean[4],
                            quantiles$mean[5]- quantiles$mean[4]))
  
  
#permutation test p values
combinations <- list(c("SR", "ST"), c("SR", "UNBIAS"), c("ST", "UNBIAS"))

perms <- 1000
outs <- list()
for (combs in combinations){
  print(combs)
  for (chrm in c("Chr_A", "Chr_X")){
    print(chrm)
    for (i in 1:perms){
      data_tmp <- filter(omega_genes, DGE_bias %in% combs & 
                       chr == chrm) %>% 
        group_by(DGE_bias) %>% 
        sample_n(min(table(.$DGE_bias)), replace = T)
       

      data_tmp$group <- sample(data_tmp$DGE_bias, replace = F)
      
      table(data_tmp$group, data_tmp$DGE_bias)
      
      all_out_data_tmp <- data_tmp %>% 
        group_by(group) %>% 
        summarise(S = sum(S, na.rm = T), 
                  N = sum(N, na.rm = T),
                  SdS = sum(SdS, na.rm = T),
                  NdN = sum(NdN, na.rm = T),
                  dS = sum(dS, na.rm = T), 
                  dN = sum(dN, na.rm = T),
                  dS2 = SdS/S,
                  dN2 = NdN/N) %>% 
        mutate(dNdS = dN/dS, 
               dNdS2 = dN2/dS2)
      
      out <- data.frame(dNdS_g1 <- all_out_data_tmp$dNdS[1], 
                        dNdS_g2 <- all_out_data_tmp$dNdS[2]) %>% 
        mutate(dif = dNdS_g1 - dNdS_g2, 
               comb = paste(combs, collapse = "_"), 
               chr = chrm)
      colnames(out) <- c("dNdS1", "dNdS2", "dif", "comb", 'chr')
      outs[[paste0(combs[1],combs[2], i, chrm)]] <- out
    }
  }
}

perm_out <- outs %>% 
  bind_rows() %>% 
  group_by(comb, chr) %>% 
  reframe(Q2.5 = quantile(dif, probs = .025), 
            Q97.5 = quantile(dif, probs = 0.975),
            dif = dif) %>% 
  unique() %>% 
  merge(lines)
  


perm_out %>% 
  mutate(comb = gsub("_", "vs", comb)) %>% 
  ggplot(aes(x = dif, color = comb)) + 
  geom_density() + 
  facet_wrap(~comb + chr, ncol = 2) +
  geom_vline(aes(xintercept = true_dif, colour = comb)) + 
  geom_vline(aes(xintercept = Q2.5, colour = comb),  linetype = 'dashed', colour = 'black') + 
  geom_vline(aes(xintercept = Q97.5, colour = comb),  linetype = 'dashed', colour = 'black')


p_values <- 
  perm_out %>% 
  group_by(comb, chr) %>% 
  reframe(more = length(which(dif > true_dif)), 
          less = length(which(dif < true_dif)),
          p = 2*(more/perms)) %>% 
  ungroup() %>% 
  mutate(group1 = str_split(comb, "_", simplify= T)[,1], 
         group2 = str_split(comb, "_", simplify = T)[,2], 
         y.position = c(0.04, 0.04, 0.04,0.04, 0.04, 0.04), 
         sig = ifelse(p < 0.01, "**", 
                      ifelse(p < 0.05, "*", 
                      "ns"))) %>% 
  mutate(fdr = p.adjust(p, method = 'fdr'),
         sig_fdr = ifelse(fdr < 0.01, "**", 
                          ifelse(fdr < 0.05, "*", 
                                 "ns"))) %>% 
  mutate(group1 = ifelse(group1 != "UNBIAS", paste0(group1, "-biased"), "Unbiased")) %>% 
  mutate(group2 = ifelse(group2 != "UNBIAS", paste0(group2, "-biased"), "Unbiased"))


p_values
boostrap_data %>% 
#  filter(chr == "Chr_A") %>% 
  mutate(class = factor(class, levels = c("SR", "ST", "UNBIAS"))) %>% 
  mutate(class = ifelse(class != "UNBIAS", paste0(class, "-biased"), "Unbiased")) %>% 
  ggplot(aes(x = values, colour = class)) +
  geom_density() + labs(x = 'difference in dNdS')
  facet_wrap(~chr)


A <- boostrap_data %>% 
  filter(chr == "Chr_A") %>% 
  mutate(class = factor(class, levels = c("SR", "ST", "UNBIAS"))) %>% 
  mutate(class = ifelse(class != "UNBIAS", paste0(class, "-biased"), "Unbiased")) %>% 
  ggplot(aes(x = class, y = values, fill = class)) + 
  geom_boxplot(notch = T, outlier.alpha = 0.2) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_classic() +
  labs(fill = "", x = "", y = 'dN/dS', title = "Autosomes")+
  stat_pvalue_manual(p_values[p_values$chr == "Chr_A"&p_values$p < 0.05,], label = 'sig') + 
  ylim(0,0.06)

X <- boostrap_data %>% 
  filter(chr == "Chr_X") %>% 
  mutate(class = factor(class, levels = c("SR", "ST", "UNBIAS"))) %>% 
  mutate(class = ifelse(class != "UNBIAS", paste0(class, "-biased"), 'Unbiased')) %>% 
  ggplot(aes(x = class, y = values, fill = class)) + 
 # geom_boxplot(notch = T, outlier.alpha = 0.2) + 
  scale_fill_brewer(palette = "Set2") + 
  theme_classic() +
  labs(fill = "", x = "", y = 'dN/dS', title = "X chromosome")+
  stat_pvalue_manual(p_values[p_values$chr == "Chr_X" &p_values$p < 0.05,], label = 'sig') + 
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(notch = T, outlier.alpha = 0.2)+ 
  ylim(0,0.06)


ggarrange(A,X, common.legend = T, align = 'h')



## ALISON IDEA 

perms <- 1000
outs <- list()
for (combs in combinations){
  print(combs)
  for (chrm in c("Chr_A", "Chr_X")){
    print(chrm)
    for (i in 1:perms){
      data_tmp <- filter(omega_genes, DGE_bias %in% combs & 
                           chr == chrm) %>% 
        group_by(DGE_bias) %>% 
        sample_n(min(table(.$DGE_bias)))
      
      
      data_tmp$group <- sample(data_tmp$DGE_bias, replace = F)
      
      table(data_tmp$group, data_tmp$DGE_bias)
      
      all_out_data_tmp <- data_tmp %>% 
        group_by(group) %>% 
        summarise(S = sum(S, na.rm = T), 
                  N = sum(N, na.rm = T),
                  SdS = sum(SdS, na.rm = T),
                  NdN = sum(NdN, na.rm = T),
                  dS = sum(dS, na.rm = T), 
                  dN = sum(dN, na.rm = T),
                  dS2 = SdS/S,
                  dN2 = NdN/N) %>% 
        mutate(dNdS = dN/dS, 
               dNdS2 = dN2/dS2)
      
      out <- data.frame(dNdS_g1 <- all_out_data_tmp$dNdS[1], 
                        dNdS_g2 <- all_out_data_tmp$dNdS[2]) %>% 
        mutate(dif = dNdS_g1 - dNdS_g2, 
               comb = paste(combs, collapse = "_"), 
               chr = chrm)
      colnames(out) <- c("dNdS1", "dNdS2", "dif", "comb", 'chr')
      outs[[paste0(combs[1],combs[2], i, chrm)]] <- out
    }
  }
}


