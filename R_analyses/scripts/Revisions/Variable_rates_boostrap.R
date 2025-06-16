library(tidyverse)
library(ape)
library(stringr)
library(plyranges) 
library(rtracklayer)
library(ggpubr)
library(boot)

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
  species_N <- sum(as.numeric(ssdata$N))
  species_S <- sum(as.numeric(ssdata$S))
  species_NdN <- sum(as.numeric(ssdata$`N*dN`))
  species_SdS <- sum(as.numeric(ssdata$`S*dS`))
  
  species_dS <- species_SdS/species_S
  species_dN <- species_NdN/species_N
  dNdS2 <- species_dN/species_dS
  tstv <- x_lines[grep("ts/tv", x_lines)] %>% strsplit(" ") %>% 
    .[[1]] %>% 
    .[length(.)] %>% as.numeric()
  
  out_data <- data.frame(species = sp, species_sci = sp_name, dS = species_dS, 
                         dN = species_dN, NdN = species_NdN, SdS = species_SdS, N = species_N, S = species_S,
                         dNdS = species_dNdS, dNdS2 <- dNdS2, orthogroup = orthogroup, 
                         tstv = tstv)
  return(out_data)
}
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
#  mutate(Significant = ifelse(logFC < -1, "SR-biased", 
#                               ifelse(logFC > 1, "ST-biased", "Unbiased"))) %>% 
  group_by(genes, chr, consensus_gene) %>% 
  summarise(ST_bias = sum(Significant == "ST-biased") > 0, 
            SR_bias = sum(Significant == "SR-biased") > 0, 
            logFC_DGE = mean(logFC)) %>% 
  mutate(DGE_bias = ifelse(ST_bias == T & SR_bias == T, "BOTH", 
                           ifelse(ST_bias == T, "ST", 
                                  ifelse(SR_bias == T, "SR", "UNBIAS")))) %>% 
  dplyr::select(-c("ST_bias", "SR_bias"))


############ COMPARING RESULTS #################

omega_genes_save <- list()
p_values_save <- list()
boostrap_data_save <- list()
comps_save <- list()
bs_summary <- list()
types <- c("broader_phy") #, "run1", "run2")
for (type in types){
  branch_folders <- list.files(paste0("data/var_rates/", type, "/"), full.names = T,  
                               pattern = '_branch_models')
  branch_files <- list.files(branch_folders, full.names = T, recursive = T)
  omegas <- lapply(branch_files, get_omega_func) %>% 
    bind_rows()
  rownames(omegas) <- NULL
  
  
  orthogroups <-read.table(paste0("data/var_rates/", type, "/N0.tsv"), sep = "\t", header = T)
  omega_genes <- merge(omegas, orthogroups, by.x = "orthogroup", by.y = 'OG', all.x = T, all.y = F) %>% 
    merge(gene_names, by.x = 'Teleopsis_dalmanni.protein_longest', by.y = 'longest_transcript') %>% 
    mutate(genes = gsub("_", "-", genes)) %>% 
    merge(dif_exp_sum, by = 'genes') %>% 
    mutate(chr2 = ifelse(chr == "Chr_X", "Chr_X", "Chr_A"), 
           chr = ifelse(chr == "Chr_X", "Chr_X", "Chr_A")) %>% 
   filter(dS <=2) %>% 
    mutate(dNdS2 = 
             (NdN/N)/(SdS/S))
  
  omega_genes_save[[type]] <- omega_genes
  #Confidence intervals 
  
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
    return(out$dNdS2)
  }
  
  
  SR_biased <- dNdS_func(omega_genes, which(omega_genes$DGE_bias == "SR"))
  ST_biased <- dNdS_func(omega_genes, which(omega_genes$DGE_bias == "ST"))
  Unbiased <-  dNdS_func(omega_genes, which(omega_genes$DGE_bias == "UNBIAS"))
 # Both <-  dNdS_func(omega_genes, which(omega_genes$DGE_bias == "BOTH"))
  
  output <- list()
  #for (Chrom in unique(omega_genes$chr2)){
    for (class in unique(omega_genes$DGE_bias)){
      print(class)
      ss <- filter(omega_genes, DGE_bias == class)# & chr2 == Chrom)
      boot_dist <- boot(data = ss, statistic = dNdS_func, R = 10000)
      boot_dist$class <- class
     # boot_dist$chr <- Chrom
      output[[paste0(class)]] <- boot_dist
      
    }
  #}
  
  
  bs_data_func <- function(boot_dist){
    out <- data.frame(values = boot_dist$t, 
               #      chr = boot_dist$chr, 
                      class = boot_dist$class)
    return(out)
  }
  
  
  boostrap_data <- lapply(output, bs_data_func) %>% 
    bind_rows()
  boostrap_data_save[[type]] <- boostrap_data
  combinations <- list(c("SR", "ST"), c("SR", "UNBIAS"), c("ST", "UNBIAS"))
  
  bs_summary[[type]] <- 
    boostrap_data %>% 
      group_by(class) %>% 
      reframe(Q2.5 = quantile(values, probs = .025), 
              Q97.5 = quantile(values, probs = 0.975), 
              mean = mean(values), 
              min = min(values), 
              max = max(values)) %>% 
      mutate(observed_mean = c(dNdS_func(omega_genes, which(omega_genes$DGE_bias == "SR")), 
                               dNdS_func(omega_genes, which(omega_genes$DGE_bias == "ST")), 
                               dNdS_func(omega_genes, which(omega_genes$DGE_bias == "UNBIAS")))) %>% 
      mutate(meanCI = paste0(class, ", mean = ", 
                          round(observed_mean, 3), 
                          " CI 95% [", 
                          round(Q2.5,3), 
                          ",", 
                          round(Q97.5,3), 
                          "]")) %>% 
      dplyr::select(meanCI)
    
  
  comps <- list()
  for (combs in combinations){
    print(combs)
    #for (chrm in c("Chr_A", "Chr_X")){
     # print(chrm)
      bs_ss <- filter(boostrap_data, class %in% combs)# & chr == chrm)
      data1 <- filter(bs_ss, class == combs[1])
      data2 <- filter(bs_ss, class == combs[2])
      difs <- outer(data1$values, data2$values, "-")
      obs_mean_diff <- mean(difs)
      out <- data.frame(Q2.5 = quantile(difs, probs = .025), 
                        Q97.5 = quantile(difs, probs = 0.975), 
                        mean = mean(difs), 
                        min = min(difs), 
                        max = max(difs), 
                     #   chr = chrm, 
                        class = paste(combs, collapse = "_"), 
                        p1 = 1-(length(difs[difs > 0])/length(difs)),
                        p2 = 1-(length(difs[difs < 0])/length(difs)), 
                        p3 = 1-(length(abs(difs)))
      ) %>% 
        mutate(p_two_tailed = 2 *min(c(p1,p2))
  )
      
      comps[[paste0(paste(combs, collapse = "_"))]] <- out
        # data1 %>% 
        #mutate(values = values - data2$values) %>% 
        #mutate(class = paste(combs, collapse = "_"))
  
    }
  #}
  
  # Print the result (a 4x4 matrix)
  
  
  comps <- comps %>% 
    bind_rows()         
  #View(comps)
  row.names(comps) <- NULL
  comps
  comps_save[[type]] <- comps
  
  p_values <- comps %>% 
    mutate(p = p_two_tailed) %>% 
    mutate(fdr = p.adjust(p, "fdr")) %>% 
    mutate(group1 = str_split(class, "_", simplify= T)[,1], 
           group2 = str_split(class, "_", simplify = T)[,2], 
           y.position = c(0.04), 
           sig = ifelse(p< 0.001, "***", ifelse(p < 0.01, "**", 
                        ifelse(p < 0.05, "*", 
                               "ns")))) %>% 
    mutate(fdr = p.adjust(p, method = 'fdr'),
           sig_fdr = ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", 
                            ifelse(fdr < 0.05, "*", 
                                   "ns")))) %>% 
    mutate(group1 = ifelse(group1 != "UNBIAS", paste0(group1, "-biased"), "Unbiased")) %>% 
    mutate(group2 = ifelse(group2 != "UNBIAS", paste0(group2, "-biased"), "Unbiased")) %>% 
   # mutate(y.position = c(0.095, 0.055, 0.085,0.045, 0.12, 0.1)) %>% 
    dplyr::rename(comp = class)
  
  
  p_values_save[[type]] <- p_values
  
  
}
  
  
boostrap_data_save$broader_phy %>% ggplot(aes(x = values, colour = class)) + geom_density()













#####################################################################################
#############################
A <- boostrap_data %>%
#  filter(chr == "Chr_A") %>% 
  mutate(class = factor(class, levels = c("SR", "ST", "UNBIAS"))) %>% 
  mutate(class = ifelse(class != "UNBIAS", paste0(class, "-biased"), "Unbiased")) %>% 
  ggplot(aes(x = class, y = values, fill = class)) + 
  geom_boxplot(notch = T, outlier.alpha = 0.2) + 
  scale_fill_manual(values = c("ST-biased" = "#66C2A5", 
                               "SR-biased" = "#FC8D62", 
                               "Unbiased" = "#8DA0CB")) + 
  theme_classic() + 
  labs(fill = "", x = "", y = 'bootstrap dN/dS')+
  stat_pvalue_manual(p_values[p_values$chr == "Chr_A"&p_values$p < 0.05,], label = 'sig_fdr') + 
  #stat_pvalue_manual(p_values[p_values$p < 0.05,], label = 'sig_fdr', y.position = 0.4) + 
 # ylim(0,0.1) + 
  theme(legend.position = 'none')

A

X <- boostrap_data %>%
  filter(chr == "Chr_X") %>% 
  mutate(class = factor(class, levels = c("SR", "ST", "UNBIAS"))) %>% 
  mutate(class = ifelse(class != "UNBIAS", paste0(class, "-biased"), "Unbiased")) %>% 
  ggplot(aes(x = class, y = (values), fill = class)) + 
  geom_boxplot(notch = T, outlier.alpha = 0.2) + 
  scale_fill_manual(values = c("ST-biased" = "#66C2A5", 
                               "SR-biased" = "#FC8D62", 
                               "Unbiased" = "#8DA0CB")) + 
  theme_classic() + 
  labs(fill = "", x = "", y = 'bootstrap dN/dS')+
  stat_pvalue_manual(p_values[p_values$chr == "Chr_X"&p_values$p < 0.05,], label = 'sig_fdr', y.position = 4) + 
  theme(legend.position = 'none')# +
 # ylim(0,0.1)
X
dNdS_plot <- ggarrange(A,X, common.legend = F, labels = c("(a)", "(b)"), label.x = -0.038)
ggsave("plots/manuscript_plots/SX_dNdS_boostrap_revisions.pdf", dNdS_plot, height = 3.5, width = 6)
system("open plots/manuscript_plots/SX_dNdS_boostrap_revisions.pdf")





