#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

install.packages(c("ape"), lib = ".")

library(ape, lib.loc = '.')
#library(ape)
#trees <- read.nexus("data/VARIABLE_RATES/TREE_17102023/output.nex")
#mcc <- maxCladeCred(trees)
#mcc$node.label <- NULL
#write.tree(mcc, "data/VARIABLE_RATES/TREE_17102023/species_tree.nex") # output in newick format


tree <- read.tree(args[1])
#tree <- read.tree("data/VARIABLE_RATES/TREE_17102023/species_tree.txt")
paml_tree <- unroot(tree)
paml_tree$edge.length <- NULL
paml_tree$node.label <- NULL


#Abreviate species names to first character of first word and first four of second word, seperated by _
for (i in 1:length(paml_tree$tip.label)){
  split <- strsplit(paml_tree$tip.label[i], "_")[[1]]
  paml_tree$tip.label[i] <- paste(substr(split[1], 1, 1), substr(split[2], 1, 3), sep = "_")
}

#paml_tree_name <- gsub(".txt", "_paml.txt", "data/VARIABLE_RATES/TREE_17102023/species_tree.txt")
paml_tree_name <- gsub(".txt", "_paml.txt", args[1])

write.tree(paml_tree, paml_tree_name)
paml_tree <- readLines(paml_tree_name)
first_row <- paste(length(tree$tip.label), "1")
paml_tree <- c(first_row, paml_tree)
writeLines(paml_tree, paml_tree_name)

