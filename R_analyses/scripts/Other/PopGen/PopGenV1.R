library(tidyverse)

depth <- read.table("indata/PopGen/sr2.thresholds.bed.gz")
colnames(depth) <- c("chr", "start", "end", "gene", "1X", "2X", '5X', "10X")
depth %>% 
  ggplot(aes(x = start, y = `5X`)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  facet_wrap(~chr, scales = "free")


Stalkie_reseq_metadata <- read.csv("indata/PopGen/pom2022_reseq_sampleinfo.csv")
Stalkie_reseq_metadata %>% 
  ggplot(aes(x = stream, fill = genotype)) + 
  geom_bar(position = 'dodge') + 
  facet_wrap(~collection.year)

Stalkie_reseq_metadata %>% 
  filter(stream == 'UBW') %>% 
  dplyr::select(sample_number, genotype) %>% 
  write.csv("indata/PopGen/UBW_genotype.csv", row.names = FALSE, quote  = F)


chrom_sizes <- read.csv("indata/PopGen/chromsizes.csv", header =F)

out <- data.frame(chr = NA, start = NA, end = NA, max = NA)
for (i in 1:3){
  seq <- seq.int(1,chrom_sizes[i,2], chrom_sizes[i,2]/25) %>% round()
  seq2 <- seq[-1] -1 %>% round()
  seq2 <- c(seq2, chrom_sizes[i,2])
  out_tmp <- data.frame(chr = chrom_sizes[i,1], start = seq, end = seq2, max = chrom_sizes[i,2])
  out <- rbind(out, out_tmp)
}
mt <- c(chrom_sizes[4,1], 1, chrom_sizes[4,2], chrom_sizes[4,2])
out <- rbind(out ,mt)
out <- out[-1,]
write.csv(out, "indata/PopGen/chrom_sizes_bins.csv", row.names = FALSE, quote = F, col.names = FALSE)

