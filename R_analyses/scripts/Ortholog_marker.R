library(tidyverse)
library(stringr)
library(readxl)
#read in a gtf file 
gtf <- read.table("indata/ref_files/OMAnnotation_output/Tdal_ST_oma.gff", sep = '\t')

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

syn_db_whole <- read.table("indata/ref_files/fb_synonym_fb_2023_06.tsv", sep = '\t', fill = T, quote = "")
colnames(syn_db_whole) <- c("primary_FBid",  "organism_abbreviation",   "current_symbol",  "current_fullname",  
              "fullname_synonym",     "symbol_synonym")

syn_db <- filter(syn_db_whole, organism_abbreviation == "Dmel") %>% #remove all rows where current_symbol ends in -RA or -RB or ]
  filter(!grepl("-RA", current_symbol)) %>%
  filter(!grepl("-RB", current_symbol)) %>%
  filter(!grepl("\\[.*\\]", current_symbol)) %>% 
  filter(current_fullname != "")

#Get synonym function

get_syn <- function(x, syn_db){
  CG <- gsub("Dmel_", "", x)
  
  #search symbol synonym column for exact match with GC
  ssmatch = grep(CG, syn_db$symbol_synonym)
  gene_name = CG
  if (length(ssmatch) == 0){
    return(data.frame(Dmel = x, CG = CG, gene = gene_name, FBgn = NA))
  } else {
    gene_pos <- syn_db %>% filter(grepl(CG, symbol_synonym))
    checks <- strsplit(gene_pos$symbol_synonym, "|", fixed = T)
    FBgn <- gene_pos$primary_FBid
    for (i in 1:length(checks)){
      if (CG %in% checks[[i]]){
        gene_name <- gene_pos$current_symbol[i]
        FBgn <- FBgn[i]
        
      }
    }
    return(data.frame(Dmel = x, CG = CG, gene = gene_name, FBgn = FBgn))
  }
}

syn_df <- genes[is.na(genes$RefGene) == F,2] %>%
  lapply(., get_syn, syn_db = syn_db) %>% 
  bind_rows()  
rownames(syn_df) <- NULL

gene_syn <- merge(genes, syn_df, by.x = "RefGene", by.y = "Dmel", all = T) %>% 
  unique()

#### READING IN ORTHOFINDER RESULTS OF DROS VS STALKIE
ortho_dir <- '/Users/peter/Documents/PhD/Projects/2022/Meiotic_drive_2022/MeioticDrive2022_Analyses/CL_analysis/orthologs_id/dros_stalkie_orth/ref_files/longest_pep'

of <- read.table(paste(ortho_dir, "//Phylogenetic_Hierarchical_Orthogroups/N0.tsv", sep = ""), sep = '\t', header = T)

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


ofk_syn <- merge(ofk, syn_db_whole[,c('primary_FBid', 'current_symbol')], by.x = "Dmel.pep_primary", by.y = "primary_FBid", all.x = T) %>% 
  dplyr::select(Tdal_ST_oma_primary, current_symbol, Dmel.pep_primary) %>% 
  unique() %>% 
  dplyr::rename(ortho_Tdal = Tdal_ST_oma_primary, ortho_Dmel = current_symbol,
                FBgnOF = Dmel.pep_primary)

final_orths <- merge(ofk_syn, gene_syn, by.x = 'ortho_Tdal', by.y = 'GeneName', all = T) %>% 
  dplyr::rename(REF_GENE_NAME = ortho_Tdal, 
         OMA_DMEL = gene,
         OMA_CG = CG,
         OMA_REFGENE = RefGene,
         OF_DMEL = ortho_Dmel, 
         FBgn_OMA = FBgn
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
  write.csv("outdata/orthologs_Jan24.csv", quote = F, row.names = F)
final_orthslwc <- final_orths %>% 
  mutate(tolower(consensus_gene))

### ELIFE MARKERS
markers_elifetmp <- readxl::read_excel("indata/markers/elife2019/elife-47138-supp1-v1.xlsx") %>% 
  mutate(Gene = tolower(Gene))
markers_elife <- final_orths %>% 
  mutate(consensus_gene = tolower(consensus_gene)) %>% 
  merge(markers_elifetmp, by.x = 'consensus_gene', by.y = 'Gene', all.y = T)

  
markers_elife_out <- markers_elife %>% 
  filter(is.na(REF_GENE_NAME) == F) %>%
  dplyr::select(Cluster, REF_GENE_NAME) %>% 
  dplyr::rename(celltype = Cluster, 
         marker = REF_GENE_NAME) %>% #replave all spaces in celltype with _
  mutate(celltype = str_replace_all(celltype, " ", "_")) %>%
  mutate(celltype = str_replace_all(celltype, ",", "_or"))

write.csv(markers_elife_out, "outdata/markers_elife.csv", quote = F, row.names = F)

######## ATLAS MARKERS ----
markers_atlas <- read_excel("indata/markers/flyatlas_dros_markers.xlsx") %>% 
  filter(Tissue == "testis")

collapse <- function(x){
  gene <- t(str_split(x[15], ",", simplify = TRUE))
  gene <- gsub(" ", "", gene, fixed = TRUE)
  gene <- c(str_split(gene, ":", simplify = TRUE))
  gene <- gene[gene != ""]
  Cluster = rep(x[2], length(gene))
  return(data.frame())
}

markers_atlaslist <- apply(markers_atlas, 1, function(x)(return(
  data.frame(
    Cluster = rep(x[2]), 
    gene = t(str_split(x[15], ",", simplify = TRUE)) %>% 
      gsub(pattern = " ", replacement =  "")) 
)))


markers_atlaslistdf <- bind_rows(markers_atlaslist) %>% 
  mutate(gene = tolower(gene))


markers_atlaslistdf$Cluster[grep('spermatoc', markers_atlaslistdf$Cluster)] <- "Spermatocytes"
markers_atlaslistdf$Cluster[grep('cyst', markers_atlaslistdf$Cluster)] <- "Cyst"
markers_atlaslistdf$Cluster[grep('epith', markers_atlaslistdf$Cluster)] <- "Epithelial cells"
markers_atlaslistdf$Cluster[grep('hub', markers_atlaslistdf$Cluster)] <- "Hub cells"
markers_atlaslistdf$Cluster[grep('spermatog', markers_atlaslistdf$Cluster)] <- "GSC, Early spermatogonia"
markers_atlaslistdf$Cluster[grep('early elongation stage spermatid', markers_atlaslistdf$Cluster)] <- "Early spermatids"
markers_atlaslistdf$Cluster[grep('early-mid elongation-stage spermatid', markers_atlaslistdf$Cluster)] <- "Early spermatids"
markers_atlaslistdf$Cluster[grep('mid-late elongation-stage spermatid', markers_atlaslistdf$Cluster)] <- "Mature spermatids"
markers_atlaslistdf$Cluster[grep('somatic', markers_atlaslistdf$Cluster)] <- "Somatic cells"

markers_atlaslistdf_out <- markers_atlaslistdf %>% 
  merge(final_orthslwc, by.x = 'gene', by.y = 'consensus_gene', all.x = T) %>% 
  dplyr::filter(is.na(REF_GENE_NAME) == F) %>%
  dplyr::select(Cluster, REF_GENE_NAME) %>% 
  dplyr::rename(celltype = Cluster, 
                marker = REF_GENE_NAME) %>% #replave all spaces in celltype with _
  mutate(celltype = str_replace_all(celltype, " ", "_")) %>%
  mutate(celltype = str_replace_all(celltype, ",", "_or"))
write.csv(markers_atlaslistdf_out, "outdata/markers_atlas.csv", quote = F, row.names = F)


markers_elife_atlas <- markers_elife_out %>% 
  merge(markers_atlaslistdf_out, all = T)
write.csv(markers_elife_atlas, "outdata/markers_elife_atlas.csv", quote = F, row.names = F)


### EBI markers --------
EBI_atlas <- read.table("indata/markers/EBI_atlas.tsv", sep = "\t", header = T)
unique(EBI_atlas$cluster)

keep_clusters <- c("adult fat cell", "epithelial cell","muscle cell", 
                   "testis sheath epithelial cell","spermatocyte", 
                   "male accessory gland main cell",'hemocyte',
                   'early elongation stage spermatid',
                   'cyst cell of testis',
                   #'germ cell',
                   'ejaculatory bulb epithelial cell',
                   'early-mid elongation-stage spermatid',
                   'mid-late elongation-stage spermatid',
                   #'germline cell',
                   'spermatid',
                   'cyst progenitor cell',
                   'male accessory gland',
                   "spermatocyte cyst cell",
                   "ejaculatory duct epithelial cell",
                   #"testis",
                   "primary spermatocyte",
                   "semen-secreting cell of the male reproductive system",
                   "seminal vesicle",
                   "male accessory gland secondary cell",
                   #"male germline cell",
                   "spermatogonium-spermatocyte transition",
                   "mature primary spermatocyte")
EBI_markers <- EBI_atlas %>% 
  filter(cluster %in% keep_clusters)
EBI_markers$cluster_simp <- EBI_markers$cluster
EBI_markers$cluster_simp[grep("spermatocyte", EBI_markers$cluster)] <- "spermatocyte"
EBI_markers$cluster_simp[grep("spermatid", EBI_markers$cluster)] <- "spermatid" #now cyst
EBI_markers$cluster_simp[grep("cyst", EBI_markers$cluster)] <- "cyst" #not epithel
EBI_markers$cluster_simp[grep("epithelial", EBI_markers$cluster)] <- "epithelial" #accessiry gland
EBI_markers$cluster_simp[grep("accessory", EBI_markers$cluster)] <- "accessory gland" 



EBI_markers_out <- EBI_markers  %>%
  merge(final_orths, by.x = 'genes', by.y = 'FBconcensus', all.x = T) %>% 
  dplyr::filter(is.na(REF_GENE_NAME) == F) %>%
  #dplyr::select(cluster_simp, REF_GENE_NAME) %>% 
  dplyr::rename(celltype = cluster_simp, 
                marker = REF_GENE_NAME) %>% #replave all spaces in celltype with _
  mutate(celltype = str_replace_all(celltype, " ", "_")) %>%
  mutate(celltype = str_replace_all(celltype, ",", "_or"))
dim(EBI_markers_out)
EBI_markers_out %>% 
  dplyr::select(cluster, marker) %>% 
  rename(celltype = cluster) %>%
  write.csv("outdata/markers_EBI.csv", quote = F, row.names = F)              

head(EBI_markers_out)


read.csv("indata/markers/cell_cycle_markers.csv") %>% 
  merge(final_orths, by.x = 'geneID', by.y = 'FBconcensus') %>% 
  dplyr::select('phase', 'geneID', 'modified', 'REF_GENE_NAME', 'consensus_gene') %>% 
  rename(gene = REF_GENE_NAME) %>% 
  write.csv('outdata/cell_cycle_markers_Jan24.csv', quote = F, row.names = F)


Mah_markers <- readxl::read_excel('indata/markers/Mahadevaraju2021.xlsx', sheet = 2, col_names = T, skip = 1) %>%
  merge(final_orths, by.x = 'FBgn', by.y = 'FBconcensus', all.x = T)

colnames(Mah_markers) <- gsub(" ", "_", colnames(Mah_markers))
Mah_markers_out <- data.frame(celltype = NULL, marker = NULL, Fb = NULL)
for (i in 1:nrow(Mah_markers)) {
  if (Mah_markers$FBgn[i] %in% final_orths$FBconcensus) {
    celltypes <- str_split(Mah_markers$Higher_Expressed_in[i], ", ", simplify = T) %>% t()
    out_data <- data.frame(celltype = celltypes, marker = Mah_markers$REF_GENE_NAME[i], Fb = Mah_markers$FBgn[i])
    Mah_markers_out <- rbind(Mah_markers_out, out_data)
  }
}

Mah_markers_out <- unique(Mah_markers_out)
rownames(Mah_markers_out) <- NULL

write.csv(Mah_markers_out, "outdata/markers_Mah.csv", quote = F, row.names = F)

          