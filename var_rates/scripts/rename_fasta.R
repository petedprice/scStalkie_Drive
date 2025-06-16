args = commandArgs(trailingOnly=TRUE)


#Read in replacement names 
names=read.table(args[1], header=F, sep="\t")
#if names do not start with ">" add ">" at start of each name

for (c in 1:ncol(names)){
  if (substr(names[,c],1,1)!=">"){
    names[,c]=paste0(">",names[,c])
  }
  #if file ends in .protein_longest remove .protein_longest
  if (substr(names[,c],nchar(names[,c])-15,nchar(names[,c]))==".protein_longest"){
    names[,c]=substr(names[,c],1,nchar(names[,c])-16)
  }
}

#read fasta file into r using base r
fasta=readLines(args[2])

#get index of which rows start with ">"
idx = grep(">", fasta)
#replace indexed lines with names from names file
fasta[idx] = unlist((names))
writeLines(fasta, args[3])


