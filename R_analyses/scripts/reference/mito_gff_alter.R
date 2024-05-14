gff <- read.table("indata/ref_files/mito/Final_stalkie_mitogenome.gff")
gff$V9 <- gsub("color=", "", gff$V9)
out_gff <- data.frame(matrix(ncol = 9))
colnames(out_gff) <- colnames(gff)
#gff$gene <- str_split(gff$V9, ";", simplify = T)[,2] %>% 
#  gsub("Name=", "", .)

for (i in 1:nrow(gff)){
  if (gff[i,3] %in% c("gene", "CDS")){
    out_gff <- rbind(out_gff, gff[i,])
  } else if(gff[i,3] == "tRNA"){
    tgff <- gff[i,]
    parent <- which(gff[,4] == tgff[,4] & gff[,5] == tgff[,5] & gff[,3] == "gene")
    pg <- gff[parent,9] %>% str_split(., ";", simplify = T) %>% .[,2] %>% 
      gsub("Name=", "", .)
    tgff[,9] <- gsub("Parent=.*", paste0("Parent=", pg), tgff[,9])
    tgff[,9] <- gsub("Parent=tRNA", "Parent=gene", tgff[,9])
    exon <- tgff
    exon[,3] <- 'exon'
    exon[,9] <- gsub("gene-", "tRNA-", exon[,9])
    out_gff <- rbind(out_gff, tgff, exon)
  } else if(gff[i,3] == "rRNA"){
    rgff <- gff[i,]
    parent <- which(gff[,4] == rgff[,4] & gff[,5] == rgff[,5] & gff[,3] == "gene")
    pg <- gff[parent,9] %>% str_split(., ";", simplify = T) %>% .[,2] %>% 
      gsub("Name=", "", .)
    rgff[,9] <- gsub("Parent=.*", paste0("Parent=", pg), rgff[,9])
    rgff[,9] <- gsub("Parent=rRNA", "Parent=gene", rgff[,9])
    exon <- rgff
    exon[,3] <- 'exon'
    exon[,9] <- gsub("gene-", "rRNA-", exon[,9])
    out_gff <- rbind(out_gff, rgff, exon)
  }
}

out_gff <- out_gff[!is.na(out_gff[,1]),]
rownames(out_gff) <- NULL
View(out_gff)

write.table(out_gff, "indata/ref_files/mito/Final_stalkie_mitogenome_altered.gff",
            quote = F, sep = "\t", row.names = F, col.names = F)
