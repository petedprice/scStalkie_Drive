#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

phylip = readLines(args[1])
#phylip = readLines("data/VARIABLE_RATES/missing_E/OG0003964_codon.nogaps_t7w15_t4w9.phy")

nospecies <- strsplit(phylip[1], " ")[[1]]
nospecies <- as.numeric(nospecies[nospecies!=""][1])
sp_indx=2
nlines=(length(phylip)-1)/nospecies
sp_indx <- c(sp_indx, ((1:(nospecies-1))*nlines)+2)
species <- paste(">", phylip[sp_indx], sep="")

#Get the seqs which are intervals between the sp_indx, save each as single sequence in seqs, the last  seq is the last sp_indx to the end of the file
seqs <- c()
for (i in 1:(length(sp_indx)-1)){
  seqs <- c(seqs, paste(phylip[(sp_indx[i]+1):(sp_indx[i+1]-1)], collapse=""))
}
seqs
#Get the last seq
seqs <- c(seqs, paste(phylip[(sp_indx[length(sp_indx)]+1):length(phylip)], collapse=""))

#Interleave teh species and seqs, so sp[1], seq[1], sp[2], seq[2] etc etc
interleaved <- c()
for (i in 1:length(seqs)){
  interleaved <- c(interleaved, species[i], seqs[i])
}

out_name <- gsub(".phy", ".fa", args[1])
writeLines(interleaved, out_name)

