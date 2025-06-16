#Marker blast 
blast_ab <- read.table("data/ref_files/new_stalkie/agat_gene_seqs/marker_blast/ST_to_Mahaj_marker.txt")
blast_ba <- read.table("data/ref_files/new_stalkie/agat_gene_seqs/marker_blast/Mahaj_to_ST_marker.txt")
colnames(blast_ab) <- c("query_a", "subject_b", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blast_ba) <- c("query_b", "subject_a", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")
best_hits_ab <- blast_ab %>%
  group_by(query_a) %>%
  arrange(evalue, desc(bitscore)) %>%
  slice(1) %>%
  ungroup()

best_hits_ba <- blast_ba %>%
  group_by(query_b) %>%
  arrange(evalue, desc(bitscore)) %>%
  slice(1) %>%
  ungroup()

Y_candidatesa <- blast_ba$subject_a #%>% 
Y_candidatesb <- blast_ab$query_a 
#  gsub("_", "-", .)

filter(ortholog_table, REF_GENE_NAME %in% Y_candidates)


p1 <- DotPlot(seurat_final, features = Y_candidates)
p2 <- DimPlot(seurat_final, label = T)
markers_plots <- ggarrange(p1, p2)
markers_plots <- annotate_figure(markers_plots, top = text_grob("marker blast", 
                                      color = "red", face = "bold", size = 14))


### Primers blast ###
blast_ab <- read.table("data/ref_files/new_stalkie/agat_gene_seqs/primer_blast/ST_to_Mahaj_primers.txt")
blast_ba <- read.table("data/ref_files/new_stalkie/agat_gene_seqs/primer_blast/Mahaj_to_ST_primers.txt")
colnames(blast_ab) <- c("query_a", "subject_b", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blast_ba) <- c("query_b", "subject_a", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")
best_hits_ab <- blast_ab %>%
  group_by(query_a) %>%
  arrange(evalue, desc(evalue)) %>%
  slice(1) %>%
  ungroup()
best_hits_ba <- blast_ba %>%
  group_by(query_b) %>%
  arrange(evalue, desc(evalue)) %>%
  slice(1) %>%
  ungroup()

Y_candidates <- best_hits_ba$subject_a %>% 
  gsub("_", "-", .)


p1 <- DotPlot(seurat_final, features = Y_candidates) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

p2 <- DimPlot(seurat_final, label = T)
primer_plots <- ggarrange(p1, p2)
primer_plots <- annotate_figure(primer_plots, top = text_grob("primer blast", 
                                                                color = "red", face = "bold", size = 14))


ggsave("plots/marker_plots.pdf", markers_plots, width = 20, height = 10)
ggsave("plots/primer_plots.pdf", primer_plots,  width = 20, height = 10)



#custom ref maker 
fasta <- readLines("data/ref_files/new_stalkie/agat_gene_seqs/Mahajan_2017_Y_linked_markers.fa")
seqs <- grep(">", fasta)
other_seqs <-(1:length(fasta))[-c(seqs)]

lengths <- data.frame(seq = NA, length = NA)
c = 0
for (i in seqs){
  c=c+1
  start = seqs[c] + 1
  end = seqs[c+1]-1
  if(i == max(seqs)){
    end = length(fasta)
  }
  seq_length <- fasta[start:end] %>% paste(collapse = "") %>% 
    strsplit("") %>% .[[1]] %>% length()
  outdata <- data.frame(seq = fasta[i], length = seq_length)
  lengths <- rbind(lengths, outdata)
}

lengths <- lengths[-is.na(lengths$seq),]
gene <- gsub(">T.dalmanni-", "", lengths$seq) %>% 
  sub("\\|.*", "", .)
lengths$gene <- gene

for (i in 1:nrow(lengths)){
  gene <- lengths$gene[i]
  gtf_line <- paste0(gene, "\tunknown\texon\t1\t", lengths$length[i], 
                     "\t.\t+\t.\tgene_id \"", 
                     gene, "\"; transcript_id \"", gene, 
                     "\"; gene_name \"", gene, "\"; gene_biotype \"protein_coding\";")
  lengths$gtf[i] <- gtf_line
  
}


new_fasta <- fasta
new_fasta[seqs] <- paste0(">", lengths$gene)

writeLines(new_fasta, "data/Y_contigs/Y_genes.fasta")
writeLines(lengths$gtf, "data/Y_contigs/Y_genes.")


load("data/RData/integrated_seurat_nf200_mtr0.20_gu0_newref_Seurat5.1_df0.08_Ycontigs.RData")
DefaultAssay(seurat_integrated) <- "RNA"
other_genes <- c("gene-6002")
seurat_integrated$treatment <- "ST"
seurat_integrated$treatment[grep("sr", seurat_integrated$sample)] <- "SR"
p1 <- DotPlot(seurat_integrated, 
        features = c(other_genes, lengths$gene), group.by = 'integrated_snn_res.0.4') + 
  coord_flip()


seurat_integrated[['yexp']] <- PercentageFeatureSet(seurat_integrated, features = intersect(rownames(seurat_integrated), lengths$gene), assay = "RNA")


#p1 <- DotPlot(seurat_integrated, features = c('PB.4560', 'yexp', 'PB.1956', 'gene-6002'), group.by = 'integrated_snn_res.0.4', split.by = 'treatment')
p2 <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.4', label = T, split.by = 'treatment')
dev.off()
ggarrange(p1,p2)


p1 <- subset(seurat_integrated, treatment == "SR") %>% 
  DotPlot(., features = lengths$gene, group.by = 'integrated_snn_res.0.4')
p2 <- subset(seurat_integrated, treatment == "ST") %>% 
  DotPlot(., features = lengths$gene, group.by = 'integrated_snn_res.0.4')
p3 <- DimPlot(seurat_integrated, group.by = 'integrated_snn_res.0.4', label = T, split.by = 'treatment')
dev.off()
ggarrange(p1,p2, p3)


