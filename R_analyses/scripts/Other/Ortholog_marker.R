library(tidyverse)
library(stringr)
library(readxl)
#read in a gtf file 
gtf <- read.table("data/ref_files/OMAnnotation_output/Tdal_ST_oma.gff", sep = '\t')

get_ref_gene_func <- function(x){
  #extract from column 9 the string starting with Dmel up until the next semicolon
  chr <- x[1][1]
  start <- x[4][1]
  end <- x[5][1]
  ref_gene <- x[9] %>% 
    strsplit(., ";", fixed = T) %>% 
      unlist() %>%  #which start with 'ref-gene' %>% 
      str_extract("ref-gene.*")  %>% #filter out NAs
      #.[!is.na(.)] %>% #remove the 'ref-gene' part
      str_replace("ref-gene=gene-", "") #remove the 'ref-gene' part
  if (sum(is.na(ref_gene) == F) == 1){
    rg <- ref_gene[is.na(ref_gene) == F]
  } else {
    rg <- NA
  }
  gene_name <- x %>% 
    strsplit(., ";", fixed = T) %>% 
    unlist() %>% 
    str_extract("Parent.*") %>% 
    #.[!is.na(.)] %>%
    str_replace("Parent=", "")
  if (sum(is.na(gene_name) == F) == 1){
    gn <- gene_name[is.na(gene_name) == F]
  } else {
    gn <- NA
  }
  return(data.frame(GeneName = gn, RefGene = rg, chr = chr, start = start, 
                    end = end))
}

mrna <- gtf %>% filter(V3 == "mRNA") %>% 
  unique()

genes <- apply(mrna, 1, get_ref_gene_func) %>% 
  bind_rows() %>% 
  unique()
rownames(genes) <- NULL

syn_db_whole <- read.table("data/ref_files/fb_synonym_fb_2023_05.tsv", sep = '\t', fill = T, quote = "")
colnames(syn_db_whole) <- c("primary_FBid",  "organism_abbreviation",   "current_symbol",  "current_fullname",  
              "fullname_synonym",     "symbol_synonym")

syn_db <- filter(syn_db_whole, organism_abbreviation == "Dmel") %>% #remove all rows where current_symbol ends in -RA or -RB or ]
  filter(grepl("FBgn", primary_FBid)) #%>%
 # filter(!grepl("-RA", current_symbol)) %>% 
#  filter(!grepl("-RB", current_symbol)) %>%
#  filter(!grepl("\\[.*\\]", current_symbol))
#  filter(current_fullname != "")

#Get synonym function

get_syn <- function(x, syn_db){
  CG <- gsub("Dmel_", "", x)
  
  #search symbol synonym column for exact match with GC
  #First check to see if CG naming is the name of the gene 
  symbol_match = which(syn_db$current_symbol == CG)
  if (length(symbol_match) == 1){
    gene_name = syn_db$current_symbol[symbol_match]
    FBgn = syn_db$primary_FBid[symbol_match]
    Synonym = syn_db$symbol_synonym[symbol_match]
    return(data.frame(Dmel = x, CG = CG, gene = gene_name, FBgn = FBgn, Synonym = Synonym))
  } else {
  
    #If it's not check the synonyms
    ssmatch = grep(CG, syn_db$symbol_synonym)
    gene_name = CG
    if (length(ssmatch) == 0){
      return(data.frame(Dmel = x, CG = CG, gene = gene_name, FBgn = NA))
    } else {
      gene_pos <- syn_db %>% filter(grepl(CG, symbol_synonym))
      checks <- strsplit(gene_pos$symbol_synonym, "|", fixed = T)
      FBgns <- gene_pos$primary_FBid
      FBgn <- NA
      for (i in 1:length(checks)){
        if (CG %in% checks[[i]]){
          gene_name <- gene_pos$current_symbol[i]
          FBgn <- FBgns[i]
          Synonym = gene_pos$symbol_synonym[i]
        }
      }
      return(data.frame(Dmel = x, CG = CG, gene = gene_name, FBgn = FBgn, Synonym = Synonym))
    }
  }
}

syn_df <- genes[is.na(genes$RefGene) == F,2] %>%
  lapply(., get_syn, syn_db = syn_db) %>% 
  bind_rows()  %>% 
  unique()
rownames(syn_df) <- NULL

gene_syn <- merge(genes, syn_df, by.x = "RefGene", by.y = "Dmel", all = T) %>% 
  unique()

#### READING IN ORTHOFINDER RESULTS OF DROS VS STALKIE
ortho_dir <- '/Users/peter/Documents/Science_Work/PhD/Projects/2022/Meiotic_drive_2022/scStalkie_Drive/CL_analysis/orthologs_id/dros_stalkie_orth/ref_files/longest_pep'

of <- read.table("data/orthology/Dmel_Tdal_Orthofinder_N0.tsv", sep = '\t', header = T)

#extact reads with one to one orthologs 
check_one_to_one <- function(x){
  out1 <- x[4] %>% as.character() %>% 
    strsplit(., ',') %>% 
    unlist()
  out2 <- x[5] %>% as.character() %>% 
    strsplit(., ',') %>% 
    unlist()
  if ((length(out1) * length(out2)) == 1){
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}


KEEP <- apply(of, 1, check_one_to_one)
ofk <- of[KEEP,]


ofk_syn <- merge(ofk, syn_db_whole[,c('primary_FBid', 'current_symbol', 'symbol_synonym')], by.x = "Dmel.pep_primary", by.y = "primary_FBid", all.x = T) %>% 
  dplyr::select(Tdal_ST_oma_primary, current_symbol, Dmel.pep_primary, symbol_synonym) %>% 
  unique() %>% 
  dplyr::rename(ortho_Tdal = Tdal_ST_oma_primary, ortho_Dmel = current_symbol,
                FBgnOF = Dmel.pep_primary, 
                OF_Synonym = symbol_synonym)

final_orths <- merge(ofk_syn, gene_syn, by.x = 'ortho_Tdal', by.y = 'GeneName', all = T) %>% 
  dplyr::rename(REF_GENE_NAME = ortho_Tdal, 
         OMA_DMEL = gene,
         OMA_CG = CG,
         OMA_REFGENE = RefGene,
         OF_DMEL = ortho_Dmel, 
         FBgn_OMA = FBgn,
         OMA_synonym = Synonym
         )


final_orths$consensus_gene <- final_orths$OMA_DMEL
ortho_extras <- is.na(final_orths$consensus_gene) &
                         ((final_orths$OF_DMEL %in% final_orths$OMA_DMEL) == F)
final_orths$consensus_gene[ortho_extras] <- final_orths$OF_DMEL[ortho_extras]

final_orths$FBconcensus <- final_orths$FBgn_OMA
FB_extras <- is.na(final_orths$FBconcensus) &
  ((final_orths$FBgnOF %in% final_orths$FBgn_OMA) == F)
final_orths$FBconcensus[FB_extras] <- final_orths$FBgnOF[FB_extras]


final_orths %>% 
  write.table("outdata/orthologs_April25.tsv", sep = "\t", quote = F, row.names = F)

