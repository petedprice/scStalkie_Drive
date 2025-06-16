#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

tree <-readLines(args[1])[2]
#tree <- readLines("data/VARIABLE_RATES/TREE_17102023/species_tree_paml.txt")[2]
#split by (, ), comma and ;
species <- strsplit(tree, "[\\(\\),;]")[[1]]
species <- species[species != ""]

rst <- readLines(args[2])
#rst <- readLines("data/VARIABLE_RATES/tree_checking/mod0/OG0004454_mod0.txt")


#### GET TIP LABELS
TREE1 <- rst[grep("TREE", rst)]
TREE2 <- strsplit(TREE1, ":")[[1]][2]
TREE3 <-strsplit(TREE2, ";")[[1]][1]
TREE4 <- strsplit(TREE3, "")[[1]]
#check if character is a number
isnum <- function(x) {
  !is.na(as.numeric(x))
}
#silence warning message
options(warn=-1)
TREETIPS<- TREE4[isnum(TREE4)]
names(species) <- as.numeric(TREETIPS)


### GET EDGE LABELS 
edges <- rst[grep("lnL", rst)+1]
edges <- strsplit(edges, " ")[[1]]
edges <- edges[edges != ""]

return_edge_label_func <- function(edge){
  #split edge by .
  out <- strsplit(edge, "\\..")[[1]]
  out <- as.numeric(out)
  return(out)
}
edgedf <- as.data.frame(do.call(rbind, lapply(edges, return_edge_label_func)))


#using the names of the tips, match up which of the edges a species sits on 
#(i.e. which edge is it the parent of)
edge_names <- c()
for (i in 1:length(edgedf[,1])){
  edge_names[i] <- paste(species[edgedf[i,2]], sep = "_")
}

edgedf$edge_name <- edge_names
while(length(c(which(edgedf$edge_name == "NA"))) > 0){
  for (edge in c(which(edgedf$edge_name == "NA"))){
    print(edge)
    parents <- which(edgedf$V1 == edgedf$V2[edge])
    species <- unique(edgedf[parents,3])
    check_NAs <- which(species == "NA")
    #remove any NA in species, return NA
    if (sum(check_NAs) > 0) {
      edgedf$edge_name[edge] <- "NA"
    } else {
      edgedf$edge_name[edge] <- paste(species, collapse = ",")
    }
  }
}

edgedf$numbering <- paste(edgedf$V1, "..", edgedf$V2, sep = "")

#write to text file with no quotes etc.
write.table(edgedf[,c(4,3)], "SWAMP_BRANCHES.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

