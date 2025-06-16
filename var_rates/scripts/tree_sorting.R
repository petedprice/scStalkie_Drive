library(ape)
library(phangorn)
trees <- read.nexus("data/VARIABLE_RATES/TREE_17102023/output.nex")
mcc <- maxCladeCred(trees)
write.tree(mcc, "data/VARIABLE_RATES/TREE_17102023/species_tree.txt") # output in newick format

mcc <- unroot(mcc)
mcc$edge.length <- NULL
mcc$node.label <- NULL


#Abreviate species names to first character of first word and first four of second word, seperated by _
for (i in 1:length(mcc$tip.label)){
  split <- strsplit(mcc$tip.label[i], "_")[[1]]
  mcc$tip.label[i] <- paste(substr(split[1], 1, 1), substr(split[2], 1, 3), sep = "_")
}

plot(mcc)
edgelabels(mcc$edge[,2])

write.tree(mcc, "data/VARIABLE_RATES/TREE_17102023/paml_species_tree.txt")
paml_tree <- readLines("data/VARIABLE_RATES/TREE_17102023/paml_species_tree.txt")
first_row <- paste(length(mcc$tip.label), "1")
paml_tree <- c(first_row, paml_tree)
writeLines(paml_tree, "data/VARIABLE_RATES/TREE_17102023/paml_species_tree.txt")

tips <- mcc$tip.label
#number each tip
names(tips) <- c(1:length(tips))

#of the reamining edges that don't touch a tip, work out which species it conects to
edges <- as.data.frame(mcc$edge)
#using the names of the tips, match up which of the edges a species sits on 
#(i.e. which edge is it the parent of)
edge_names <- c()
for (i in 1:length(edges[,1])){
  edge_names[i] <- paste(tips[edges[i,2]], sep = "_")
}

edges$edge_name <- edge_names
while(length(c(which(edges$edge_name == "NA"))) > 0){
  for (edge in c(which(edges$edge_name == "NA"))){
    print(edge)
    parents <- which(edges$V1 == edges$V2[edge])
    species <- unique(edges[parents,3])
    check_NAs <- which(species == "NA")
    #remove any NA in species, return NA
    if (sum(check_NAs) > 0) {
      edges$edge_name[edge] <- "NA"
    } else {
      edges$edge_name[edge] <- paste(species, collapse = ",")
    }
  }
}

edges$numbering <- paste(edges$V1, "..", edges$V2, sep = "")

#write to text file with no quotes etc.
write.table(edges[,c(4,3)], "data/VARIABLE_RATES/TREE_17102023/edges.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')

