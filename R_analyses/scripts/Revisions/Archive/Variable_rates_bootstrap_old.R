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



omega_genes <- omega_genes %>% 
  mutate(chr2 = ifelse(chr == "Chr_X", "Chr_X", "Chr_A"), 
         chr = ifelse(chr == "Chr_X", "Chr_X", "Chr_A"))
#Confidence intervals 

combinations <- list(c("SR", "ST"), c("SR", "UNBIAS"), c("ST", "UNBIAS"))

dNdS_func <- function(df){
  out <- df %>% 
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
for (Chrom in unique(omega_genes$chr2)){
  print(Chrom)
  for (class in unique(omega_genes$DGE_bias)){
    print(class)
    for (B in 1:100){
      ss <- filter(omega_genes, chr2 == Chrom & DGE_bias == class)
      g1 <- ss %>% 
        sample_n(10)
      dNdS <- dNdS_func(g1)
      out <- data.frame(dNdS = dNdS, class = class, chr = Chrom)
      
      
      
      output[[paste0(class, Chrom, B)]] <- out
    }  
  }
}

output <- output %>% 
  bind_rows()


comps <- list()
for (combs in combinations){
  print(combs)
  for (chrm in c("Chr_A", "Chr_X")){
    print(chrm)
    bs_ss <- filter(output, chr == chrm & class %in% combs)
    data1 <- filter(bs_ss, class == combs[1])
    data2 <- filter(bs_ss, class == combs[2])
    difs <- outer(data1$dNdS, data2$dNdS, "-") %>% c()
    out <- data.frame(Q2.5 = quantile(difs, probs = .025), 
                      Q97.5 = quantile(difs, probs = 0.975), 
                      mean = mean(difs), 
                      min = min(difs), 
                      max = max(difs), 
                      chr = chrm, 
                      class = paste(combs, collapse = "_"))
    
    
    comps[[paste0(chrm, paste(combs, collapse = "_"))]] <- out
    # data1 %>% 
    #mutate(values = values - data2$values) %>% 
    #mutate(class = paste(combs, collapse = "_"))
    
  }
}
comps <- comps %>% 
  bind_rows()
row.names(comps) <- NULL
comps
