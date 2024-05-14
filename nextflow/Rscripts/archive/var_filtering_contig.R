#!/usr/bin/env Rscript
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

snp_info <- read.table(args[1], comment.char = "") %>% 
  as.data.frame()

colnames(snp_info) <- c("Read-Name", "Barcode",  "Flag",  "MAPQ", "CHROM",  
                        "READ-POS0",  "READ-BASE",  "READ-QUAL",  "REF-POS1",
                        "REF-BASE",  "CIGAR-OP")

snp_info$`REF-POS1` <- suppressWarnings(as.numeric(snp_info$`REF-POS1`))
snp_info <- snp_info[is.na(snp_info$`REF-POS1`) == F,]
snp_info$`READ-QUAL` <- sapply(snp_info$`READ-QUAL`, utf8ToInt) 
snp_info$`READ-BASE` <- toupper(snp_info$`READ-BASE`)
snp_info$`REF-BASE` <- toupper(snp_info$`REF-BASE`)

snp_info <- snp_info %>% 
  filter(`READ-QUAL` > 62)

snp_summarise <- snp_info %>% 
  group_by(Barcode, `REF-POS1`) %>% 
  reframe(refn = length(which(`READ-BASE` == `REF-BASE`)), 
          altn = length(which(`READ-BASE` != `REF-BASE`)), 
          CHROM = `CHROM`[1], 
          ref = `REF-BASE`[1], 
          alt = `READ-BASE`[1],
          alt = alt[1]) %>% 
  rename(pos = `REF-POS1`)


write.table(snp_summarise, file = gzfile(paste(args[2], args[3], "snp_summarise.txt.gz", sep = "_")), row.names = F, 
	col.names = F, quote = F)
