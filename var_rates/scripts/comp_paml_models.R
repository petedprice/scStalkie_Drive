library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
files=list.files(args[1], pattern = "paml", full.names = T)
folders=list.dirs(args[1], full.names = T, recursive = F)

#files=list.files("data/VARIABLE_RATES/delete", full.names = T)
#folders=list.dirs("data/VARIABLE_RATES/delete", full.names = T, recursive = F)

files = files[(files %in% folders)==F]
nmods=args[2]


#files=list.files("data/VARIABLE_RATES/paml_output/", full.names = T)
#for each file return the likelihoods, the line starts with lnL
#and the model name, the line starts with Model
if (length(files) > 0){
  model_data <- data.frame()
  duds <- c()
  for (i in 1:length(files)){
    print(i)
    file <- files[i]
    file_data <- readLines(file)
    lnL_row <- file_data[grep("lnL",file_data)]
    #get np value, value after np: and before close bracket
    np <- gsub(".*np: ","",lnL_row)
    np <- gsub("\\).*","",np)
    #get rid of everything before ):
    lnL <- gsub(".*\\): ","",lnL_row)
    #split this by space
    lnL <- unlist(strsplit(lnL," "))
    #get rid of empty elements
    lnL <- lnL[lnL != ""]
    #keep every other from first
    lnL <- lnL[c(T,F)]
    if (length(lnL) != nmods){
      duds <- c(duds, file)
      next("Number of models does not match number of lnL values")
    }
    
    model_line <- file_data[grep("Model ",file_data)]
    model <- gsub("NSsites ","",model_line)
    model <- gsub(":.*","",model)
    model <- gsub(" ","_",model)
    
    orthogroup <- gsub(".*OG","",file)
    orthogroup <- paste("OG", gsub("_.*","",orthogroup), sep = "")
    
    tmp_md <- as.data.frame(
      matrix(nrow = 1, ncol = (length(lnL)*2), data = c(as.numeric(np),as.numeric(lnL)))
    )
    colnames(tmp_md) <- c(paste(model, "np", sep="_"), paste(model, "lnL", sep="_"))
    tmp_md$orthogroup <- orthogroup
    model_data <- rbind(model_data, tmp_md)
  }
  model_data
  model_data$ln_dif <- 2*(model_data$Model_2_lnL - model_data$Model_1_lnL)
  model_data$pval <- unlist(lapply(model_data$ln_dif, pchisq, lower.tail = F, df = 2))
  #sort table with lowest p value first
  model_data <- model_data[order(model_data$pval),]
  
  write.table(model_data, paste0("m1a2a", args[3]), row.names = F, col.names = T, quote = F)
  write.table(duds, paste("duds_m1a2a_", args[3], sep = ""), row.names = F, col.names = F, quote = F)
}  


ln_val_func <- function(file){
  file_data <- readLines(file)
  lnL_row <- file_data[grep("lnL",file_data)]
  #get np value, value after np: and before close bracket
  np <- gsub(".*np: ","",lnL_row)
  np <- gsub("\\).*","",np)
  #get rid of everything before ):
  lnL <- gsub(".*\\): ","",lnL_row)
  #split this by space
  lnL <- unlist(strsplit(lnL," "))
  #get rid of empty elements
  lnL <- lnL[lnL != ""]
  #keep every other from first
  lnL <- lnL[c(T,F)]
  branch1 <- strsplit(file, "/")[[1]]
  branch <- strsplit(branch1[length(branch1)], "_")[[1]][1]
  if (length(grep("null", branch1)) ==1){
    model <- "null_branch"
  } else {
    model <- "branch"
  }
  
  return(list(np = np, lnL = as.numeric(lnL), branch = branch, model = model))
}


comp_model_func <- function(branch, model_data){
  sm_data = model_data[model_data$branch == branch,]
  sm_data$ln_dif <- 2*(sm_data$lnL[1] - sm_data$lnL[2])
  sm_data$pval <- pchisq(sm_data$ln_dif, lower.tail = F, df = 1)
  return(sm_data)
}


if (length(folders) > 0){
  model_data <- data.frame()
  duds <- c()
  for (i in 1:length(folders)){
    print(i)
    folder <- folders[i]
    ffiles <- list.files(folder, pattern = "paml", full.names = T)
    nmods <- length(ffiles)
    tmp_md <- bind_rows(lapply(ffiles, ln_val_func))
    tmp_md2 <- bind_rows(lapply(
      unique(tmp_md$branch), comp_model_func, tmp_md))
    if (length(tmp_md2$lnL) != nmods){
      duds <- c(duds, file)
      next("Number of models does not match number of lnL values")
    }
    
    orthogroup <- gsub(".*OG","",ffiles[1])
    orthogroup <- paste("OG", gsub("_.*","",orthogroup), sep = "")
    tmp_md2$orthogroup <- orthogroup
    model_data <- rbind(model_data, tmp_md2)
  }

  #sort table with lowest p value first
  model_data <- model_data[order(model_data$pval),]
  
  write.table(model_data, paste0("branch_", args[3]), row.names = F, col.names = T, quote = F)
  write.table(duds, paste("duds_branch_", args[3], sep = ""), row.names = F, col.names = F, quote = F)
  
}

