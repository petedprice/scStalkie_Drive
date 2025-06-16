#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

install.packages(c("ape"), lib = ".")
library(ape, lib.loc = '.')

trees <- list.files(args[1], full.names = T, pattern = 'rooted')
dir.create("paml_branch_trees")

tree_branch_function <- function(x){
  tree <- read.tree(x)
  paml_tree <- unroot(tree)
  paml_tree$edge.length <- NULL
  
  #Abreviate species names to first character of first word and first four of second word, seperated by _
  for (i in 1:length(paml_tree$tip.label)){
    name=paml_tree$tip.label[i]
    split1 <- strsplit(name, "#")[[1]]
    split2 <- strsplit(split1, "_")[[1]]
    newname <- paste(substr(split2[1], 1, 1), substr(split2[2], 1, 3), sep = "_")
    if (length(split1) > 1 ){
      newname <- paste(newname, split1[2], sep = "#")
    }
    print(newname)
    paml_tree$tip.label[i]  <- newname
  }
  outfile1 = strsplit(x, "/")[[1]]
  outfile2 = paste0("paml_branch_trees/", outfile1[length(outfile1)])
  paml_tree_name <- gsub("rooted_branch.txt", "paml_branch.txt", outfile2)
  
  write.tree(paml_tree, paml_tree_name)
  paml_tree <- readLines(paml_tree_name)
  first_row <- paste(length(tree$tip.label), "1")
  paml_tree <- c(first_row, paml_tree)
  writeLines(paml_tree, paml_tree_name)
}

# make directory 
lapply(trees, tree_branch_function)

