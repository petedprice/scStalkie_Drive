library(stringr)
library(tidyverse)

mummer_data <- readLines("indata/ref_files/ST_SR_X.mums")
md2 <- str_split(mummer_data, " ")
x <- md2[2]

mum_data_get <- function(x){
  x2 <- x %>% 
    unlist() %>% 
    as.numeric() 
  x2 <- x2[which(!is.na(x2))]
  if (length(x2) == 3){
    return(x2)
  } else {
    return(c(0,0,0))
  }
}

plot_data <- lapply(md2, mum_data_get) %>% #glue rows together
  unlist() %>% #unlist
  matrix(ncol = 3, byrow = TRUE) %>% #make matrix
  as.data.frame() %>% #make data frame
  setNames(c("start", "end", "length"))  %>% #name columns 
  filter(start != 0 & end != 0)

dim(plot_data)
plot <-  plot_data %>% 
  ggplot(aes(x = start, y = end)) + geom_point(alpha = 0.005)


ggsave(filename = 'plots/mummer.png', plot)



#### comparing GTFs ----
get_ref_gene_func <- function(x){
  #extract from column 9 the string starting with Dmel up until the next semicolon
  chr <- x[1]
  start <- x[4]
  end <- x[5]
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

gtf_st <- read.table("indata/ref_files/OMAnnotation_output/Tdal_ST_oma.gff", sep = '\t') %>% 
  filter(V3 == "mRNA") %>% 
  unique() %>% 
  apply(1, get_ref_gene_func) %>% 
  bind_rows()
rownames(gtf_st) <- NULL

gtf_sr <- read.table("indata/ref_files/OMAnnotation_output/Tdal_SR_oma.gff", sep = '\t') %>% 
  filter(V3 == "mRNA") %>% 
  unique() %>% 
  apply(1, get_ref_gene_func) %>% 
  bind_rows()

rownames(gtf_sr) <- NULL

X_st <- gtf_st %>% 
  filter(chr == "Chr_X")
X_sr <- gtf_sr %>%
  filter(chr == "PGA_scaffold_3__41_contigs__length_92227466")
X_sr$chr <- "Chr_X"

X_st <- filter(X_st, RefGene %in% X_sr$RefGene & is.na(RefGene) == F) %>% 
  unique()
X_sr <- filter(X_sr, RefGene %in% X_st$RefGene & is.na(RefGene) == F) %>% 
  unique()

X_comb <- merge(X_st, X_sr, by = "RefGene", suffixes = c("_st", "_sr"))

gtf_cords <- X_comb %>% 
  ggplot(aes(x = as.numeric(start_st), y = as.numeric(start_sr))) + geom_point()

ggsave(gtf_cords, filename = "plots/gtf_cords.png")

#ID inversion boundries 
X_comb <- X_comb[order(X_comb$start_st),]
head(X_comb)
gtf_st$inversions <- 0
X_comb$sign <- -1

inversions <- c()

for (i in 2:(nrow(X_comb))){
  sr = X_comb$start_sr %>% 
    as.numeric()
  
  drct_sr <- sign(sr[i] - sr[i-1])
  X_comb$sign[i] <- drct_sr
  
}  

for (i in 2:(nrow(X_comb)-5)){
  signs <- X_comb$sign[i:(i+5)]
  msign <- median(signs)
  X_comb$inversions[i] <- msign
}  


X_comb$inversions <- 0
for (i in 10:(nrow(X_comb)-10)){
  b <- median(X_comb$sign[i:(i-10)])
  a <- median(X_comb$sign[i:(i+10)])
  
  if (a != b & a %in% c(1,-1) & b %in% c(1,-1)){
    X_comb$inversions[i:nrow(X_comb)] <- X_comb$inversions[i:nrow(X_comb)] +1
    print(i)
  }
}



X_comb %>% 
  ggplot(aes(x = as.numeric(start_st), y = as.numeric(start_sr), colour = as.factor(inversions))) +
  geom_point() + 
  facet_wrap(~inversions)
  #scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

