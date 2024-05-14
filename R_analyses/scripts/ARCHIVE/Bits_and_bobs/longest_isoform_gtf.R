library(dplyr)
library(ape)
library(tidyr)
library(stringr)
library(future.apply)
library(future)
library(optparse)
option_list = list(
  make_option(c("-g", "--gff"), type="character", default=".", 
              help="gff file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$gff)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

gff <- read.gff(opt$gff)
gff <- filter(gff, type %in% c("gene", "CDS"))

gene_name_func <- function(x){
    input=str_split(x, ";")
  i=which(startsWith(input[[1]], "gene="))
  return(gsub("gene=", "", input[[1]][i]))
}

gff$gene_name <- unlist(lapply(gff$attributes, gene_name_func))


#### FUNCTION FOR MAKING PSEUDO CDS FILE 
unique_cds_func <- function(g, gff){
  ss_gff <- filter(gff, gene_name == g)
  if (nrow(ss_gff) == 1){
    return(NULL)
  } else {
    cds <- apply(ss_gff[2:nrow(ss_gff),], 1, function(x)(return(x[4]:x[5]))) %>% 
      unlist() %>% 
      unique() %>% 
      sort()
    
    swaps1 <- c(unlist(lapply(2:length(cds), function(x)(if((cds[x]-cds[x-1]) > 1)(return(x))))))
    starts <- cds[c(1,swaps1)]
    
    swaps2 <- unlist(lapply(1:(length(cds)-1), function(x)(if((cds[x]-cds[x+1]) < -1)(return(x)))))
    ends <- cds[c(swaps2, length(cds))]
    #swaps2+ ss_gff$start[1] %>% sort()
    
    
    
    out_gff = ss_gff[1:(length(ends)+1),]
    out_gff[2:nrow(out_gff),] <- out_gff[2,]
    out_gff[2:nrow(out_gff),c(4,5)] <- matrix(c(starts, ends), ncol = 2)
    return(out_gff)
  }
}



### RUNNING AND CHECKING 
print("running get longest pseudo cds region")
gff_out <- lapply(unique(gff$gene_name)[1:100], unique_cds_func, gff = gff) %>% bind_rows()
out_name = gsub(".gff", "_longest.gff", opt$gff)
print("saving file")
gff_out[,c(6,8)] = "."
phase <- lapply(gff_out$start[gff_out$type == "CDS"], function(x)(return((gff$phase[gff$start == x & gff$type == "CDS"][1])))) %>%
  unlist() %>% 
  as.numeric_version() %>% 
  unlist()
gff_out$phase[gff_out$type == "CDS"]<- phase

write.table(gff_out[,-ncol(gff_out)], out_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
       