#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

phy_file = args[1]
penalty=as.numeric(args[2])
#phy_file = "data/VARIABLE_RATES/NoNs/OG0003964_codon.nogaps_t7w15_t4w9_NoNs.phy"
phy = readLines(phy_file)
#get species which is the lines with names and not the A,C,T,G,N lines
seq_info_line <- phy[1]
nspecies <- strsplit(seq_info_line, " ")[[1]]
nspecies <- as.numeric(nspecies[nspecies != ""][1])
gaps <- (length(phy)-1)/nspecies
sp_indx <- c(2,(c(1:(nspecies-1)) * gaps)+2)
species <- phy[sp_indx]
N_indx <- c()

seqs <- list()


for (i in 1:length(species)){
  si <- sp_indx[i]
  seq1 <- paste(phy[(si+1):(si+gaps-1)], collapse = "")
  seq2 <- strsplit(seq1, "")[[1]]
  #split seq by N but keep Ns
  chunks <- unlist(strsplit(seq1, "(?<=N)(?=[ATCG])|(?<=[ATCG])(?=N)", perl = TRUE))
  #check if all chunks divisible by 3
  if (sum(nchar(chunks) %%3) > 0){
    stop("Ns or seqs not divisible by 3")
  }
  #see which letters in seq are N
  N_tmp <- which(seq2 == "N")
  N_indx <- c(N_indx, N_tmp)
  seqs[[i]] <- seq2
}
species <- paste(species, "  ", sep = "")
remove_ns <- function(x, N_indx){
  x[N_indx] <- ""
  x <- paste(x, collapse = "")
  return(x)
}

seqs_NoN <-lapply(seqs, remove_ns, N_indx = N_indx)

### checks 
lengths <- lapply(seqs_NoN, nchar)
if (length(unique(lengths)) > 1){
  stop("seqs not same length")
}

if (lengths[1] < penalty){
	stop("seqs less than penalty length")
}

NuNs <- lapply(seqs_NoN, function(x) sum(strsplit(x, "")[[1]] == "N"))
#if any Ns then stop 
if (sum(unlist(NuNs)) > 0){
  stop("Ns still in seqs")
}

out_file <- strsplit(phy[1], " ")[[1]]

out_file <- out_file[out_file != ""]
out_file[2] <- lengths[1]
out_file <- paste(out_file, collapse = " ")
for (i in 1:length(species)){
  out_file <- c(out_file, species[i])
  #split seqs into groups of 60
  seqs_NoN_tmp <- strsplit(seqs_NoN[[i]], "")
  seqs_NoN_tmp <- unlist(seqs_NoN_tmp)
  seqs_NoN_tmp <- split(seqs_NoN_tmp, ceiling(seq_along(seqs_NoN_tmp)/60))
  seqs_NoN_tmp <- unlist(lapply(seqs_NoN_tmp, paste, collapse = ""))
  out_file <- c(out_file, seqs_NoN_tmp)
}

  
#writeLines(out_file, "data/VARIABLE_RATES/test_NoNs.phy")


#change output name by putting a NoN before the .phy
out_name <- gsub(".phy", "_NoNs.phy", phy_file)
writeLines(out_file, out_name)














