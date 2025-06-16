#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

allignment <- readLines(args[1])
penalty <- as.numeric(args[2])
#which lines start with >
species_idx <- which(substr(allignment, 1, 1) == ">")

species_idx <- c(species_idx, length(allignment)+1)
species <- allignment[species_idx][1:(length(species_idx)-1)]

sequences <- allignment[-species_idx]


seqs <- character()
gaps <- list()
for (i in 1:length(species)){
  seq <- paste(allignment[(species_idx[i]+1):(species_idx[i+1]-1)], collapse = "")
  seqs[i] <- seq
  #check which character index is a gap
  gaps[[i]] <- which(strsplit(seq, "")[[1]] == "-")
}

#Remove intersecting gap indexes from each sequence 
unique_gaps <- unlist(unique(gaps))
if (length(unique_gaps) > 0){
  for (s in 1:length(seqs)){
    seqs[(s)] <- paste(strsplit(seqs[(s)], "")[[1]][-unique_gaps], collapse = "")
  }
}

#Check if any sequence ends in triplet of TAA, TAG, or TGA and return true or false vector 
ends_in_stop <- c()
for (i in 1:length(seqs)){
  ends_in_stop[i] <- grepl("TAA|TAG|TGA", substr(seqs[i], nchar(seqs[i])-2, nchar(seqs[i])))
}

#If any sequence ends in stop, remove last three nucltodies from all sequences (not the species names though)
if (sum(ends_in_stop) > 0){
  for (i in 1:length(seqs)){
    seqs[i] <- substr(seqs[i], 1, nchar(seqs[i])-3)
  }
}

#Check all sequences same length and divisable by three. Return true or false vector 
lengths <- c()
for (i in 1:length(seqs)){
  lengths[i] <- nchar(seqs[i])
}

divis3 <- lengths%%3 == 0
if (sum(divis3) != length(divis3)){
  stop("Sequences not all divisible by 3")
}

#check all sequences same length
if (length(unique(lengths)) > 1){
  stop("Sequences not all same length")
} 

#interweave seqs and species
final_fasta <- c()
for (i in 1:length(species)){
  final_fasta <- c(final_fasta, species[i], seqs[i])
}

print(paste("final sequence length is ", nchar(final_fasta[2])))
print(paste("penalty limit is ", penalty))

if (nchar(final_fasta[2]) < penalty){
  stop("sequence too short")
}

writeLines(final_fasta, args[3])

